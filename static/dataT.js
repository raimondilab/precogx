function attempt(assay, assays, path_to_json_output) {
  var dataT = initTable(path_to_json_output);
  var string = [];
  //alert(assays.split(',')[2]);
  for (let i = 0; i < assays.split(',').length; i++) {
    if (document.getElementById(assays.split(',')[i]).checked == false) {
      string.push(assays.split(',')[i]);
      //alert();
    }
  }
  string = string.join('|');
  var regex = '';
  if (string != []) {
    var regex = '^(?:(?!'+string+').)*$';
  }
  dataT.columns(1).search(regex, true, false).draw();
}

function makeDatatable(path_to_json_output, path_to_fasta, uniq_id, gpcr_list) {
  $(document).ready(function() {
    //var dataT = initTable();
    initTable(path_to_json_output);
    var assays = ['IUPHAR', 'LogRAi', 'Emax', 'WT'];
    //alert(path_to_json_output);
    var checkB = document.createElement("DIV");
    checkB.setAttribute("class", "checkboxRight");
    for (let i = 0; i < assays.length; i++) {
      //alert(i);
      var check = document.createElement("DIV");
      check.setAttribute("class", "form-check form-check-inline");
      var inP = document.createElement("INPUT");
      inP.setAttribute("class", "form-check-input");
      inP.setAttribute("type", "checkbox");
      inP.setAttribute("id", assays[i]);
      inP.setAttribute("value", "option1");
      inP.setAttribute("checked", "");
      //inP.setAttribute("onchange", "alert();");
      inP.setAttribute("onclick", "attempt(\'"+assays[i]+"\',\'"+assays+"\', \'"+path_to_json_output+"\');");
      //inP.onclick = attempt();
      var label = document.createElement("LABEL");
      label.setAttribute("class", "form-check-label");
      label.setAttribute("for", assays[i]);
      label.innerHTML = assays[i];
      check.appendChild(inP);
      check.appendChild(label);
      checkB.appendChild(check);
    }
    $("div.toolbar").html(checkB);
    //$("#exportFiles").html('Hello');
    //var exportFiles = document.getElementById('exportFiles');
    //alert(exportFiles.innerHTML);

    //firstrow();
    $('#example tbody').on( 'click', 'td', function () {
      var dataT = initTable(path_to_json_output);
      //alert(dataT.cell( this ).data());
      //alert(dataT.row(this).index());
      $('#example tbody td').css('backgroundColor', 'white').data('gpcr', gpcr).data('gprotein', gprotein);
      $('#example tbody td').attr('id', 'white').css( "border", "None" );
      //alert ($('#example tbody td').html());
      var colIndex = dataT.column(this).index();
      //alert(colIndex);
      var cell = dataT.cell(this).node();
      $(cell).css('backgroundColor', 'darkgrey').css( "border", "3px solid black" ).attr('id', 'selected');
      var rowIndex = dataT.row(this).index();
      header = ['GPCR', 'VAR', 'GNAS', 'GNAL', 'GNAI1', 'GNAI2', 'GNAI3', 'GoA', 'GoB', 'GNAZ', 'GNA11', 'GNA14', 'GNA15', 'GNAQ', 'GNA12', 'GNA13', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2'];
      //alert (header[index]);
      //var gpcrs = ['DRD4'];
      var gpcrs = gpcr_list;
      var gprotein = header[colIndex];
      var gpcr = gpcrs[rowIndex];
      var variant = gpcr_list[rowIndex];
      //alert(gpcr);
      $(cell).data('gpcr', gpcr);
      $(cell).data('gprotein', gprotein);
      var slider1_value = document.getElementById('slider1_value').innerHTML;
      //alert(document.getElementById('slider1_value').innerHTML);
      makeSequence(gpcr, path_to_fasta, gprotein, slider1_value, uniq_id);
      makeStructure(gpcr, gprotein, slider1_value, uniq_id);
      makeHeatmap(slider1_value, gpcr, gprotein);
      makePCA(uniq_id, 'Shedding', 'GPCRome', gpcr, gprotein);
    });
  });
}

function firstrow(path_to_json_output) {
  //alert('first row');
  var dataT = initTable(path_to_json_output);
  var rowIndex = dataT.row(0).index();
  var colIndex = 14;
  var cell = dataT.cell(rowIndex,colIndex).node();
  $(cell).css('backgroundColor', 'darkgrey').css( "border", "3px solid black" ).attr('id', 'selected');
  var gpcr_variant = dataT.cell(0,0).data() + '_';
  gpcr_variant += dataT.cell(0,1).data();
  //alert(d);
  //$(cell).data('gpcr', dataT.cell(0,0).data());
  $(cell).data('gpcr', gpcr_variant);
  $(cell).data('gprotein', 'GNA12');
}

function initTable(path_to_json_output) {
  return $('#example').DataTable({
    //"ajax": "static/OL820/out.json",
    "ajax": path_to_json_output,
    "retrieve": true,
    "initComplete": function(  ) {
                    //alert('hello');
                    firstrow(path_to_json_output);
                    },
    "drawCallback": function( settings ) {
                    //alert('hi');
                  },
    "scrollX": true,
    "scrollY": "250px",
    "scrollCollapse": true,
    "paging": false,
    "bSort": false,
    //"dom": 'Bfrtip',
    "dom": '<"toolbar">Bfrtip',
    "buttons": [
            {
                extend: 'spacer',
                style: 'bar',
                text: '<div>G-protein groups'
            },
            {
                extend: 'colvisGroup',
                text: 'Gs',
                show: [ 0, 1, 2, 3],
                hide: [ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ]
            },
            {
                extend: 'colvisGroup',
                text: 'Gi/Go',
                show: [ 0, 1, 4, 5, 6, 7, 8, 9],
                hide: [ 2, 3, 10, 11, 12, 13, 14, 15, 16, 17, 18 ]
            },
            {
                extend: 'colvisGroup',
                text: 'Gq/G11',
                show: [ 0, 1, 10, 11, 12, 13],
                hide: [ 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 18 ]
            },
            {
                extend: 'colvisGroup',
                text: 'G12/G13',
                show: [ 0, 1, 14, 15],
                hide: [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18 ]
            },
            {
                extend: 'colvisGroup',
                text: 'B-arrs',
                show: [ 0, 1, 16, 17, 18],
                hide: [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ]
            },
            {
                extend: 'colvisGroup',
                text: 'Show all',
                show: [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ],
                hide: [ ]
            },
            {
                extend: 'spacer',
                style: 'bar',
                text: '<div id="exportFiles"> Export files'
            },
            'copy',
            'print',
            'csv',
            'excel',
            'spacer',
            {
                extend: 'spacer',
                style: 'bar',
                text: '</div>',
            },
             {
                extend: 'pdfHtml5',
                orientation: 'landscape',
                pageSize: 'LEGAL'
            }
        ],
    "fixedColumns":   {
            left: 2
        },

    columns: [
        { title: "GPCR" },
        { title: "VAR" },
        { title: "GNAS" },
        { title: "GNAL" },
        { title: "GNAI1" },
        { title: "GNAI2" },
        { title: "GNAI3" },
        { title: "GoA" },
        { title: "GoB" },
        { title: "GNAZ" },
        { title: "GNA11" },
        { title: "GNA14" },
        { title: "GNA15" },
        { title: "GNAQ" },
        { title: "GNA12" },
        { title: "GNA13" },
        { title: "Barr1-GRK2" },
        { title: "Barr2" },
        { title: "Barr2-GRK2" }
      ]
  });
}
