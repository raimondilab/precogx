
function makeDatatable(path_to_json_output, path_to_fasta, uniq_id, gpcr_list) {
  $(document).ready(function() {
    //var dataT = initTable();
    initTable(path_to_json_output);
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
      header = ['GPCR', 'GNAL', 'GNAI2', 'GNA12', 'GoA', 'GNA15', 'GNAI1', 'GNAI3', 'GNAS', 'GNAQ', 'GNA14', 'GNA11', 'GNAO1', 'GNAZ', 'GNA13', 'GoB'];
      //alert (header[index]);
      //var gpcrs = ['DRD4'];
      var gpcrs = gpcr_list;
      var gprotein = header[colIndex];
      var gpcr = gpcrs[rowIndex];
      //alert(rowIndex);
      //alert(dataT.row(0).data());
      $(cell).data('gpcr', gpcr);
      $(cell).data('gprotein', gprotein);
      var slider1_value = document.getElementById('slider1_value').innerHTML;
      //alert(document.getElementById('slider1_value').innerHTML);
      makeSequence(gpcr, path_to_fasta, gprotein, slider1_value, uniq_id);
      makeStructure(gpcr, gprotein, slider1_value, uniq_id);
      makeHeatmap(slider1_value, gpcr, gprotein);
      makePCA(uniq_id, 'Shedding', gpcr, gprotein);
    });
  });
}

function firstrow(path_to_json_output) {
  //alert('first row');
  var dataT = initTable(path_to_json_output);
  var rowIndex = dataT.row(0).index();
  var colIndex = 3;
  var cell = dataT.cell(rowIndex,colIndex).node();
  $(cell).css('backgroundColor', 'darkgrey').css( "border", "3px solid black" ).attr('id', 'selected');
  $(cell).data('gpcr', dataT.cell(0,0).data());
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

    columns: [
        { title: "GPCR" },
        { title: "GNAL" },
        { title: "GNAI2" },
        { title: "GNA12" },
        { title: "GoA" },
        { title: "GNA15" },
        { title: "GNAI1" },
        { title: "GNAI3" },
        { title: "GNAS" },
        { title: "GNAQ" },
        { title: "GNA14" },
        { title: "GNA11" },
        { title: "GNAO1" },
        { title: "GNAZ" },
        { title: "GNA13" },
        { title: "GoB" }
      ]
  });
}
