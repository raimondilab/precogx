$( document ).ready(function() {

var selected_gpcr = '';


$(function(){
    $("#modalCoupling").load("/static/html/modals.html");

});


<!-- Load the first input and GNA12 pair-->
var path_to_json_output = document.getElementById('path_to_json_output').value;
var path_to_fasta = document.getElementById('path_to_fasta').value;
var uniq_id = document.getElementById('uniq_id').value;
var gpcr_list = document.getElementById('gpcr_list').value;
var first_gprotein = document.getElementById('first_gprotein').value;
var first_gprotein_index = document.getElementById('first_gprotein_index').value;
var first_entry = document.getElementById('first_entry').value;

console.log(path_to_fasta);

makeDatatable(path_to_json_output, path_to_fasta, uniq_id, gpcr_list, first_gprotein, first_gprotein_index);
makeStructure(first_entry, first_gprotein, 0.0, 0.0, uniq_id);
makeSequence(first_entry, path_to_fasta, first_gprotein, 0.0, 0.0, uniq_id, 'Sequence');
makeHeatmap(0.0, 0.0, first_entry, first_gprotein, 'Sequence');
makeAttentionmap(uniq_id, first_entry, first_gprotein, 'Sequence');
setDisplayMenu(path_to_fasta, 0.0, 0.0, uniq_id, first_entry, first_gprotein, 'Sequence');
makePCA(uniq_id, 'TGF', '33', first_entry, first_gprotein);



//Update slider2 if changed
var slider2A_value = document.getElementById('slider2A_value');
slider2A_value.innerHTML = $('#slider2A').val();


$(document).on('input', '#slider2A', function() {
     $('#slider2A_value').html($(this).val() );
      var slider2A_value = document.getElementById('slider2A_value').innerHTML;
      var slider2B_value = document.getElementById('slider2B_value').innerHTML;
      document.getElementById('slider1A_value').innerHTML = document.getElementById('slider2A_value').innerHTML;
      var slider1A = document.getElementById('slider1A');
      slider1A.value = slider2A_value;
      var selected_gpcr = $('#selected').data('gpcr');
      var selected_gprotein = $('#selected').data('gprotein');
      var displayName = document.getElementById('displayName').innerHTML;
      var uniq_id = document.getElementById('uniq_id').value;
      var path_to_fasta = document.getElementById('path_to_fasta').value;

      makeStructure(selected_gpcr, selected_gprotein, slider2A_value, slider2B_value,uniq_id);
      makeHeatmap(slider2A_value, slider2B_value, selected_gpcr, selected_gprotein, displayName);
      makeAttentionmap(uniq_id, selected_gpcr, selected_gprotein, displayName);
      setDisplayMenu(path_to_fasta, slider2A_value, slider2B_value, uniq_id, selected_gpcr, selected_gprotein, displayName);
      makeSequence(selected_gpcr, path_to_fasta, selected_gprotein, slider2A_value, slider2B_value, uniq_id, displayName);
      makePCA(uniq_id, 'TGF', '33', selected_gpcr, selected_gprotein);
});





//End
});




function loadIntro() {
            introJs().start();
}



