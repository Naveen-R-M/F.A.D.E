/*
 * Reporting Module
 * Generates comprehensive pipeline reports
 */

process generateReport {
    tag "final_report"
    label 'process_low'
    publishDir "${params.output_dir}/07_report", mode: 'copy'
    
    container 'python:3.9'
    
    input:
    path target_info
    path structure  
    path docking_results
    path optimized_leads
    
    output:
    path "final_report.html", emit: html_report
    path "final_report.json", emit: json_report
    path "summary.txt", emit: summary
    
    script:
    """
    # Get current date
    CURRENT_DATE=\$(date)
    CURRENT_ISO=\$(date -Iseconds)
    
    # Create simple HTML report
    cat > final_report.html << HTMLEOF
    <!DOCTYPE html>
    <html>
    <head>
        <title>F.A.D.E Pipeline Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1 { color: #2c3e50; }
            .info { background: #f0f0f0; padding: 10px; margin: 10px 0; }
        </style>
    </head>
    <body>
        <h1>F.A.D.E Drug Discovery Pipeline Report</h1>
        <div class="info">
            <h2>Query</h2>
            <p>${params.query}</p>
        </div>
        <div class="info">
            <h2>Results</h2>
            <p>Pipeline completed successfully</p>
            <p>Date: \$CURRENT_DATE</p>
            <p>Target information: ${target_info}</p>
            <p>Docking results: ${docking_results}</p>
        </div>
    </body>
    </html>
    HTMLEOF
    
    # Create JSON report
    cat > final_report.json << JSONEOF
    {
        "status": "success",
        "query": "${params.query}",
        "date": "\$CURRENT_ISO",
        "files": {
            "target_info": "${target_info}",
            "structure": "${structure}",
            "docking_results": "${docking_results}",
            "optimized_leads": "${optimized_leads}"
        }
    }
    JSONEOF
    
    # Create text summary
    cat > summary.txt << SUMMARYEOF
    F.A.D.E Pipeline Summary
    ========================
    Query: ${params.query}
    Date: \$CURRENT_DATE
    Status: Complete
    
    Input Files:
    - Target: ${target_info}
    - Structure: ${structure}
    - Docking: ${docking_results}
    - Leads: ${optimized_leads}
    
    Pipeline completed successfully!
    SUMMARYEOF
    """
    
    stub:
    """
    echo "<html><body><h1>F.A.D.E Report</h1></body></html>" > final_report.html
    echo '{"status": "success"}' > final_report.json
    echo "F.A.D.E Pipeline Summary" > summary.txt
    """
}

workflow REPORTING {
    take:
    target_info
    structure
    docking_results
    optimized_leads
    
    main:
    generateReport(
        target_info,
        structure,
        docking_results,
        optimized_leads
    )
    
    emit:
    html_report = generateReport.out.html_report
    json_report = generateReport.out.json_report
    summary = generateReport.out.summary
}
