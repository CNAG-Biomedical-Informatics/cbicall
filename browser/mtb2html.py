#!/usr/bin/env python3
"""
mtb2html.py - Generate mtDNA HTML viewer (BFF-like) for MToolBox output

Usage:
  mtb2json.py -i mit_prioritized_variants.txt -f json4html > mtDNA.json
  mtb2html.py --id peter --job-id 1234 --json mtDNA.json --out mtDNA.html
"""

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate mtDNA HTML viewer for mtb2json.py json4html output"
    )
    p.add_argument("--id", required=True,
                   help="Project ID to display in the HTML header")
    p.add_argument("--job-id", required=True,
                   help="Job ID to display in the HTML header")
    p.add_argument("--json", default="mtDNA.json",
                   help="JSON file from mtb2json.py -f json4html")
    p.add_argument("--out", "-o", default="mtDNA.html",
                   help="Output HTML file name")
    return p.parse_args()


def main():
    args = parse_args()

    project_id = args.id
    job_id = args.job_id
    json_file = args.json
    html_out = args.out

    # Hardcoded extra files under ../01_mtoolbox/
    base_dir = "../01_mtoolbox"
    vcf_file = base_dir + "/VCF_file.vcf"
    report_file = base_dir + "/mit_prioritized_variants.txt"
    haplo_file = base_dir + "/mt_classification_best_results.csv"
    json_raw_file = base_dir + "/mit.raw.json"

    # Download buttons HTML
    download_buttons = "\n      ".join([
        '<a class="btn pull-right" href="{rep}"><i class="icon-download"></i> REPORT</a>'.format(rep=report_file),
        '<a class="btn pull-right" href="{hap}"><i class="icon-download"></i> HAPLOg</a>'.format(hap=haplo_file),
        '<a class="btn pull-right" href="{vcf}"><i class="icon-download"></i> VCF</a>'.format(vcf=vcf_file),
        '<a class="btn pull-right" href="./{json}"><i class="icon-download"></i> mtDNA JSON</a>'.format(json=json_raw_file),
    ])

    try:
        out = open(html_out, "w", encoding="utf-8")
    except OSError as e:
        sys.exit("Cannot write {}: {}".format(html_out, e))

    disease_threshold = 0.80  # used only for coloring

    # Use simple placeholders that we will .replace later:
    # __JSON_FILE__, __PROJECT_ID__, __JOB_ID__, __DOWNLOAD_BUTTONS__, __DISEASE_THRESHOLD__
    html = """<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>BFF Browser - mtDNA</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <link rel="icon" href="assets/img/favicon.ico" type="image/x-icon" />
    <link rel="stylesheet" href="assets/css/bootstrap.css">
    <link rel="stylesheet" href="assets/css/bootstrap-responsive.css">
    <link rel="stylesheet" href="assets/css/main.css">
    <link rel="stylesheet" href="assets/jsD/media/css/jquery.dataTables.css">
    <link rel="stylesheet" href="assets/jsD/media/css/dataTables.colReorder.css">
    <link rel="stylesheet" href="assets/jsD/media/css/dataTables.colVis.css">
    <link rel="stylesheet" href="assets/jsD/media/css/dataTables.tableTools.css">

    <script src="assets/js/jquery.min.js"></script>
    <script src="assets/js/bootstrap.min.js"></script>
    <script src="assets/jsD/media/js/jquery.dataTables.min.js"></script>
    <script src="assets/jsD/media/js/dataTables.colReorder.js"></script>
    <script src="assets/jsD/media/js/dataTables.colVis.js"></script>
    <script src="assets/jsD/media/js/dataTables.tableTools.js"></script>

    <script type="text/javascript">
    $(document).ready(function() {

      var diseaseThreshold = __DISEASE_THRESHOLD__;

      // FILTER: keep variants with ANY supporting evidence
      $.fn.dataTableExt.afnFiltering.push(
        function(oSettings, aData, iDataIndex) {

          if (oSettings.nTable.getAttribute('id') !== 'table-panel-mtdna') {
            return true;
          }
          if (!$('#toggle-candidates').prop('checked')) {
            return true;
          }

          // Columns:
          // 12 Mitomap
          // 15 ClinVar
          // 16 OMIM
          // 17 dbSNP
          var mitomap = (aData[12] || "").replace(/<[^>]*>/g, '').trim();
          var clinvar = (aData[15] || "").replace(/<[^>]*>/g, '').trim();
          var omim    = (aData[16] || "").replace(/<[^>]*>/g, '').trim();
          var dbsnp   = (aData[17] || "").replace(/<[^>]*>/g, '').trim();

          if (mitomap || clinvar || omim || dbsnp)
            return true;

          return false;
        }
      );

      // DATATABLE
      var table = $('#table-panel-mtdna').dataTable({
        "ajax": {
          "url": "__JSON_FILE__",
          "dataSrc": "data"
        },
        "bDeferRender": true,
        "stateSave": true,
        "order": [[ 8, "desc" ]],
        "search": { "regex": true },

        "language": {
          "sSearch": '<span class="icon-search"></span>',
          "lengthMenu": "Show _MENU_ variants",
          "sInfo": "Showing _START_ to _END_ of _TOTAL_ variants",
          "sInfoFiltered": " (filtered from _MAX_ variants)"
        },

        "aoColumnDefs": [
          { "aTargets": [9,11,18,13,14,20,21], "bVisible": false },

          {
            "aTargets": [10], /* Disease Score */
            "fnCreatedCell": function(td, cellData) {
              var val = parseFloat(String(cellData).replace(/<[^>]*>/g,''));
              if (!isNaN(val) && val >= diseaseThreshold)
                $(td).css('background-color','#f9d2d2');
            }
          },

          {
            "aTargets": [12], /* Mitomap */
            "fnCreatedCell": function(td, cellData) {
              var txt = String(cellData).replace(/<[^>]*>/g,'').trim();
              if (txt.length > 0)
                $(td).css('background-color','#f9d2d2');
            }
          }
        ],

        "dom": 'CRT<"clear">lfrtip',
        "colVis": {
          "showAll": "Show all",
          "showNone": "Show none"
        },
        "tableTools": {
          "aButtons": [
            {
              "sExtends": "print",
              "sButtonText": '<span class="icon-print"></span>'
            }
          ]
        },

        "orderCellsTop": true
      });


      // PER-COLUMN FILTERS
      $('#table-panel-mtdna thead tr:eq(1) th').each(function(i) {
        if (i === 1 || i === 10 || i === 15 || i === 17) {
          $(this).html(
            '<input type="text" placeholder="Filter" ' +
            'style="width:100%; box-sizing:border-box; font-size:10px;" />'
          );
        }
      });

      $('#table-panel-mtdna thead tr:eq(1) th input').on('keyup change', function() {
        var index = $(this).parent().index();
        table.fnFilter(this.value, index);
      });

      $('#toggle-candidates').on('change', function() {
        table.fnDraw();
      });

    });
    </script>

  </head>
  <body class="dt-example">

    <div class="navbar navbar-inverse navbar-fixed-top">
      <div class="navbar-inner">
        <div class="container">
          <a class="brand" href="#">BFF Browser - mtDNA</a>
        </div>
      </div>
    </div>

    <div class="container">
      __DOWNLOAD_BUTTONS__

      <h4>Project &#9658; __PROJECT_ID__</h4>
      <h3>Job ID &#9658; __JOB_ID__ &#9658; mtDNA_prioritized_variants</h3>

      <label class="checkbox">
        <input type="checkbox" id="toggle-candidates">
        Show only candidate variants (Mitomap OR ClinVar OR OMIM OR dbSNP)
      </label>

      <ul class="nav nav-tabs">
        <li class="active"><a href="#tab-panel-mtdna" data-toggle="tab">mtDNA variants</a></li>
      </ul>

      <div class="tab-content">
        <div class="tab-pane fade in active" id="tab-panel-mtdna">
          <table id="table-panel-mtdna" class="display table table-hover table-condensed">
            <thead>
              <tr>
                <th>Sample</th><th>Locus</th><th>Variant_Allele</th><th>Ref</th>
                <th>Alt</th><th>Aa_Change</th><th>GT</th><th>Depth</th>
                <th>Heterop_Frac</th><th>tRNA_Annotation</th><th>Disease_Score</th>
                <th>RNA_predictions</th><th>Mitomap_Associated_Disease(s)</th>
                <th>Mitomap_Homoplasmy</th><th>Mitomap_Heteroplasmy</th><th>ClinVar</th>
                <th>OMIM_link</th><th>dbSNP_ID</th><th>Mamit-tRNA_link</th>
                <th>AC/AN_1000G</th><th>1000G_Homoplasmy</th><th>1000G_Heteroplasmy</th>
              </tr>
              <tr>
                <th></th><th></th><th></th><th></th><th></th>
                <th></th><th></th><th></th><th></th><th></th>
                <th></th><th></th><th></th><th></th><th></th>
                <th></th><th></th><th></th><th></th><th></th>
                <th></th><th></th>
              </tr>
            </thead>
          </table>
        </div>
      </div>

      <footer><p>&copy; 2021–2025 CNAG | Barcelona</p></footer>
    </div>
  </body>
</html>
"""

    # Now substitute placeholders safely
    html = html.replace("__JSON_FILE__", json_file)
    html = html.replace("__PROJECT_ID__", project_id)
    html = html.replace("__JOB_ID__", job_id)
    html = html.replace("__DOWNLOAD_BUTTONS__", download_buttons)
    html = html.replace("__DISEASE_THRESHOLD__", "{:.2f}".format(disease_threshold))

    out.write(html)
    out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
