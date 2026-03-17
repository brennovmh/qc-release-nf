nextflow.enable.dsl=2


def samplesheetToChannel(csv_file) {
    Channel
        .fromPath(csv_file)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_fastqc.zip"), path("${sample_id}_R2_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_R1_fastqc.html"), path("${sample_id}_R2_fastqc.html"), emit: html

    script:
    """
    ln -s ${read1} ${sample_id}_R1.fastq.gz
    ln -s ${read2} ${sample_id}_R2.fastq.gz

    fastqc \
      --threads ${params.fastqc_threads} \
      ${sample_id}_R1.fastq.gz \
      ${sample_id}_R2.fastq.gz
    """
}

process FASTP_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.fastp.json"), path("${sample_id}.fastp.html"), emit: reports

    script:
    """
    fastp \
      --in1 ${read1} \
      --in2 ${read2} \
      --out1 ${sample_id}.tmp_R1.fastq.gz \
      --out2 ${sample_id}.tmp_R2.fastq.gz \
      --thread ${params.fastp_threads} \
      --json ${sample_id}.fastp.json \
      --html ${sample_id}.fastp.html

    rm -f ${sample_id}.tmp_R1.fastq.gz ${sample_id}.tmp_R2.fastq.gz
    """
}

process QC_GATE {
    tag "qc_gate"
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path fastp_jsons
    path fastqc_zips

    output:
    path "qc_summary.tsv", emit: tsv
    path "qc_summary.csv", emit: csv

    script:
    """
    python3 << 'PY'
    import json, zipfile, re, csv, glob
    from pathlib import Path

    min_reads = int(${params.min_reads})
    min_q30_rate = float(${params.min_q30_rate})
    max_adapter_percent = float(${params.max_adapter_percent})
    max_failed_fastqc_modules = int(${params.max_failed_fastqc_modules})

    fastp_files = sorted(glob.glob("*.fastp.json"))
    fastqc_files = sorted(glob.glob("*_fastqc.zip"))

    fastqc_by_sample = {}

    for z in fastqc_files:
        name = Path(z).name
        sample = re.sub(r'_R[12]_fastqc.zip\$', '', name)
        fastqc_by_sample.setdefault(sample, []).append(z)

    rows = []

    for fp in fastp_files:
        sample = Path(fp).name.replace(".fastp.json", "")
        with open(fp) as f:
            data = json.load(f)

        summary = data.get("summary", {})
        before = summary.get("before_filtering", {})
        filtering = summary.get("after_filtering", {})

        total_reads = before.get("total_reads", 0)
        q30_rate = filtering.get("q30_rate", 0)

        adapter_bases = data.get("adapter_cutting", {}).get("adapter_trimmed_bases", 0)
        total_bases = before.get("total_bases", 0)
        adapter_percent = (adapter_bases / total_bases * 100) if total_bases else 0

        failed_modules = 0
        warn_modules = 0

        for fqzip in fastqc_by_sample.get(sample, []):
            with zipfile.ZipFile(fqzip) as zf:
                summary_file = [x for x in zf.namelist() if x.endswith("/summary.txt")][0]
                content = zf.read(summary_file).decode()

            for line in content.strip().splitlines():
                status, module, filename = line.split("\\t")
                if status == "FAIL":
                    failed_modules += 1
                elif status == "WARN":
                    warn_modules += 1

        reasons = []

        if total_reads < min_reads:
            reasons.append(f"low_reads<{min_reads}")
        if q30_rate < min_q30_rate:
            reasons.append(f"low_q30<{min_q30_rate}")
        if adapter_percent > max_adapter_percent:
            reasons.append(f"high_adapter>{max_adapter_percent}%")
        if failed_modules > max_failed_fastqc_modules:
            reasons.append(f"too_many_fastqc_fail>{max_failed_fastqc_modules}")

        if failed_modules > max_failed_fastqc_modules or q30_rate < (min_q30_rate - 0.10):
            status = "FAIL"
        elif reasons or warn_modules > 7:
            status = "REVIEW"
        else:
            status = "PASS"

        rows.append({
            "sample": sample,
            "status": status,
            "total_reads": total_reads,
            "q30_rate": round(q30_rate, 4),
            "adapter_percent": round(adapter_percent, 4),
            "fastqc_fail_modules": failed_modules,
            "fastqc_warn_modules": warn_modules,
            "reasons": ";".join(reasons) if reasons else "OK"
        })

    rows = sorted(rows, key=lambda x: x["sample"])

    with open("qc_summary.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\\t")
        writer.writeheader()
        writer.writerows(rows)

    with open("qc_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    PY
    """
}

process QC_VISUAL_REPORT {
    tag "qc_visual_report"
    publishDir "${params.outdir}/report", mode: 'copy'
    container 'python:3.11-slim'

    input:
    path qc_summary_tsv

    output:
    path "run_qc_report.html", emit: html

    script:
    """
    python3 << 'PY'
    from pathlib import Path
    import csv
    import html

    infile = Path("${qc_summary_tsv}")
    rows = []
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter='\\t')
        for row in reader:
            rows.append(row)

    def badge(status):
        if status == "PASS":
            return '<span class="badge pass">PASS</span>'
        elif status == "REVIEW":
            return '<span class="badge review">REVIEW</span>'
        return '<span class="badge fail">FAIL</span>'

    def comment(row):
        status = row["status"]
        if status == "PASS":
            return "Apta para análise secundária."
        if status == "REVIEW":
            if row["reasons"] != "OK":
                return f"Requer revisão: {html.escape(row['reasons'])}."
            return "Requer revisão manual."
        return "Não apta para liberação automática."

    total = len(rows)
    n_pass = sum(r["status"] == "PASS" for r in rows)
    n_review = sum(r["status"] == "REVIEW" for r in rows)
    n_fail = sum(r["status"] == "FAIL" for r in rows)

    html_text = f'''
    <!DOCTYPE html>
    <html lang="pt-br">
    <head>
      <meta charset="utf-8">
      <title>Relatório de QC da corrida</title>
      <style>
        @page {{
          size: A4 landscape;
          margin: 12mm;
        }}

        body {{
          font-family: Arial, sans-serif;
          margin: 0;
          color: #222;
          font-size: 11px;
          line-height: 1.25;
        }}

        h1 {{
          margin: 0 0 4px 0;
          font-size: 20px;
        }}

        .subtitle {{
          color: #666;
          margin-bottom: 12px;
          font-size: 11px;
        }}

        .summary {{
          display: flex;
          gap: 8px;
          margin-bottom: 12px;
          flex-wrap: wrap;
        }}

        .card {{
          border: 1px solid #ddd;
          border-radius: 8px;
          padding: 8px 10px;
          min-width: 100px;
          background: #fafafa;
        }}

        .card .value {{
          font-size: 18px;
          font-weight: bold;
          margin-top: 2px;
        }}

        table {{
          width: 100%;
          border-collapse: collapse;
          table-layout: fixed;
          font-size: 10px;
        }}

        th, td {{
          border: 1px solid #ddd;
          padding: 6px 7px;
          text-align: left;
          vertical-align: top;
          word-wrap: break-word;
          overflow-wrap: break-word;
        }}

        th {{
          background: #f3f4f6;
          font-size: 10px;
        }}

        .badge {{
          display: inline-block;
          padding: 3px 8px;
          border-radius: 999px;
          font-weight: bold;
          font-size: 10px;
        }}

        .pass {{
          background: #d1fae5;
          color: #065f46;
        }}

        .review {{
          background: #fef3c7;
          color: #92400e;
        }}

        .fail {{
          background: #fee2e2;
          color: #991b1b;
        }}

        .small {{
          color: #666;
          font-size: 9px;
          margin-top: 10px;
        }}
      </style>
    </head>
    <body>
      <h1>Relatório de controle de qualidade da corrida</h1>
      <div class="subtitle">Resumo automatizado a partir do arquivo qc_summary.tsv</div>

      <div class="summary">
        <div class="card"><div>Total de amostras</div><div class="value">{total}</div></div>
        <div class="card"><div>PASS</div><div class="value">{n_pass}</div></div>
        <div class="card"><div>REVIEW</div><div class="value">{n_review}</div></div>
        <div class="card"><div>FAIL</div><div class="value">{n_fail}</div></div>
      </div>

      <table>
        <thead>
          <tr>
            <th>Amostra</th>
            <th>Status</th>
            <th>Total de reads</th>
            <th>Q30</th>
            <th>Adaptador (%)</th>
            <th>FastQC FAIL</th>
            <th>FastQC WARN</th>
            <th>Reasons</th>
            <th>Comentário</th>
          </tr>
        </thead>
        <tbody>
    '''

    for r in rows:
        html_text += f'''
          <tr>
            <td>{html.escape(r["sample"])}</td>
            <td>{badge(r["status"])}</td>
            <td>{html.escape(r["total_reads"])}</td>
            <td>{html.escape(r["q30_rate"])}</td>
            <td>{html.escape(r["adapter_percent"])}</td>
            <td>{html.escape(r["fastqc_fail_modules"])}</td>
            <td>{html.escape(r["fastqc_warn_modules"])}</td>
            <td>{html.escape(r["reasons"])}</td>
            <td>{comment(r)}</td>
          </tr>
        '''

    html_text += '''
        </tbody>
      </table>

      <p class="small">Este relatório é um resumo visual auxiliar. A decisão final pode considerar também o MultiQC e a revisão técnica da corrida.</p>
    </body>
    </html>
    '''

    Path("run_qc_report.html").write_text(html_text, encoding="utf-8")
    PY
    """
}
process QC_REPORT_PDF {
    tag "qc_report_pdf"
    publishDir "${params.outdir}/report", mode: 'copy'
    container 'qc-weasyprint:1.0'

    input:
    path report_html

    output:
    path "run_qc_report.pdf", emit: pdf

    script:
    """
    python3 << 'PY'
    from weasyprint import HTML
    HTML("run_qc_report.html").write_pdf("run_qc_report.pdf")
    PY
    """
}

process MULTIQC {
    tag "multiqc"
    container 'multiqc/multiqc:latest'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_files
    path fastp_files
    path qc_summary_file

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    mkdir -p multiqc_input/fastqc
    mkdir -p multiqc_input/fastp
    mkdir -p multiqc_input/summary

    cp ${fastqc_files.join(' ')} multiqc_input/fastqc/
    cp ${fastp_files.join(' ')} multiqc_input/fastp/
    cp ${qc_summary_file} multiqc_input/summary/

    multiqc multiqc_input \
      --force \
      --outdir . \
      --filename multiqc_report.html \
      -m fastqc \
      -m fastp
    """
}
workflow {
    samples_ch = samplesheetToChannel(params.input)

    fastqc_out = FASTQC(samples_ch)
    fastp_out  = FASTP_QC(samples_ch)

    fastqc_zip_list = fastqc_out.zip
        .flatMap { sample_id, zip1, zip2 -> [zip1, zip2] }
        .collect()

    fastp_json_list = fastp_out.reports
        .map { sample_id, json_file, html_file -> json_file }
        .collect()

    fastp_all_reports = fastp_out.reports
        .flatMap { sample_id, json_file, html_file -> [json_file, html_file] }
        .collect()

qc_gate = QC_GATE(fastp_json_list, fastqc_zip_list)

qc_visual = QC_VISUAL_REPORT(qc_gate.tsv)

qc_pdf    = QC_REPORT_PDF(qc_visual.html)

MULTIQC(fastqc_zip_list, fastp_all_reports, qc_gate.tsv)
}
