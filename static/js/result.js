$( document ).ready(function() {

var selected_gpcr = '';


$(function(){
    $("#modalCoupling").load("/static/html/modals.html");

});


//NGL
var stage = new NGL.Stage( "viewport" );
stage.setParameters( { backgroundColor: "white"} );

<!-- Load the first input and GNA12 pair-->
var path_to_json_output = document.getElementById('path_to_json_output').value.replace('"', "").replace('"', '');
var path_to_fasta = document.getElementById('path_to_fasta').value.replace('"', "").replace('"', '');
var uniq_id = document.getElementById('uniq_id').value.replace('"', "").replace('"', '');
var gpcr_list = document.getElementById('gpcr_list').value.replace('"', "").replace('"', '');
var first_gprotein = document.getElementById('first_gprotein').value.replace('"', "").replace('"', '');
var first_gprotein_index = document.getElementById('first_gprotein_index').value.replace('"', "").replace('"', '');
var first_entry = document.getElementById('first_entry').value.replace('"', "").replace('"', '');

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
      var selected_gpcr = $('#selected').data('gpcr').replace('"', "").replace('"', '');
      var selected_gprotein = $('#selected').data('gprotein').replace('"', "").replace('"', '');
      var displayName = document.getElementById('displayName').innerHTML;
      var uniq_id = document.getElementById('uniq_id').value.replace('"', "").replace('"', '');
      var path_to_fasta = document.getElementById('path_to_fasta').value.replace('"', "").replace('"', '');

      makeStructure(selected_gpcr, selected_gprotein, slider2A_value, slider2B_value,uniq_id);
      makeHeatmap(slider2A_value, slider2B_value, selected_gpcr, selected_gprotein, displayName);
      makeAttentionmap(uniq_id, selected_gpcr, selected_gprotein, displayName);
      setDisplayMenu(path_to_fasta, slider2A_value, slider2B_value, uniq_id, selected_gpcr, selected_gprotein, displayName);
      makeSequence(selected_gpcr, path_to_fasta, selected_gprotein, slider2A_value, slider2B_value, uniq_id, displayName);
      makePCA(uniq_id, 'TGF', '33', selected_gpcr, selected_gprotein);
});


 //Change all sliders if any of the slider button us clicked
$("#minus1A").click(function(event) {
    zoomA("out", "#slider1A", "#slider2A");
 });

$("#plus1A").click(function(event) {
    zoomA("in", "#slider1A", "#slider2A");
});

$("#minus1B").click(function(event) {
     zoomB("out", "#slider1B", "#slider2B");
});

$("#plus1B").click(function(event) {
   zoomB("in", "#slider1B", "#slider2B");
});

$("#minus2A").click(function(event) {
    zoomA("out", "#slider2A", "#slider1A");
});

$("#plus2A").click(function(event) {
  zoomA("in", "#slider2A", "#slider1A");
});

$("#minus2B").click(function(event) {
    zoomB("out", "#slider2B",  "#slider1B");
});

$("#plus2B").click(function(event) {
    zoomB("in", "#slider2B", "#slider1B");
});

//End
});



// Load intro panels
function loadIntro() {
    introJs().start();
}


// Function to reload structure panel if user selects different 3D strcture AFsource and insert in the pdblist ID
function click3Dsource(element){
    var selected_gpcr = $('#selected').data('gpcr');
    var selected_gprotein = $('#selected').data('gprotein');
    var slider2A_value = document.getElementById('slider2A_value').innerHTML;
    var slider2B_value = document.getElementById('slider2B_value').innerHTML;
    var uniq_id = document.getElementById('uniq_id').value.replace('"', "").replace('"', '');
    makeStructure(selected_gpcr, selected_gprotein, slider2A_value, slider2B_value, uniq_id);
}

//Function to update slider A i.e. log-odds score
function zoomA(direction, slider, otherSlider1, otherSlider2) {
    var slider = $(slider);
    var step = 0.1;
    var currentSliderValue = Number(slider.val());


    if (direction === "out") {
           var newStepValue = currentSliderValue - step;
    } else {
           var newStepValue = currentSliderValue + step;
    }

    newStepValue = Number(newStepValue.toFixed(1));


    slider.val(newStepValue).change();

    var sliderValue = $("#"+slider.attr("id")+'_value');
    sliderValue.html(newStepValue);

    var otherSlider1 = $(otherSlider1);
    otherSlider1.val(newStepValue).change();
    var otherSlider1Value = $("#"+otherSlider1.attr("id")+'_value');
    otherSlider1Value.html(newStepValue);

    var otherSlider2 = $(otherSlider2);
    otherSlider2.val(newStepValue).change();
    var otherSlider2Value = $("#"+otherSlider2.attr("id")+'_value');
    otherSlider2Value.html(newStepValue);

    var slider_value = slider.val();
    var selected_gpcr = $('#selected').data('gpcr').replace('"', "").replace('"', '');
    var selected_gprotein = $('#selected').data('gprotein').replace('"', "").replace('"', '');
    var slider1B_value = document.getElementById('slider1B_value').innerHTML;
    var displayName = document.getElementById('displayName').innerHTML;
    var uniq_id = document.getElementById('uniq_id').value.replace('"', "").replace('"', '');
    var path_to_fasta = document.getElementById('path_to_fasta').value.replace('"', "").replace('"', '');

    makeStructure(selected_gpcr, selected_gprotein, slider_value, slider1B_value, uniq_id);
    makeHeatmap(slider_value, slider1B_value, selected_gpcr, selected_gprotein, displayName);
    makeAttentionmap(uniq_id, selected_gpcr, selected_gprotein, displayName)
    setDisplayMenu(path_to_fasta, slider_value, slider1B_value, uniq_id, selected_gpcr, selected_gprotein, displayName);
    makeSequence(selected_gpcr, path_to_fasta, selected_gprotein, slider_value, slider1B_value, uniq_id, displayName);
    makePCA(uniq_id, 'TGF', '33', selected_gpcr, selected_gprotein);
}

//Function to update slider B i.e. distance
function zoomB(direction, slider, otherSlider1, otherSlider2) {
    var slider = $(slider);
    var step = 1;
    var currentSliderValue = Number(slider.val());


    if (direction === "out") {
         var newStepValue = currentSliderValue - step;
    } else {
         var newStepValue = currentSliderValue + step;
    }

    newStepValue = Number(newStepValue.toFixed(1));


    slider.val(newStepValue).change();

    var sliderValue = $("#"+slider.attr("id")+'_value');
    sliderValue.html(newStepValue);

    var otherSlider1 = $(otherSlider1);
    otherSlider1.val(newStepValue).change();
    var otherSlider1Value = $("#"+otherSlider1.attr("id")+'_value');
    otherSlider1Value.html(newStepValue);

    var otherSlider2 = $(otherSlider2);
    otherSlider2.val(newStepValue).change();
    var otherSlider2Value = $("#"+otherSlider2.attr("id")+'_value');
    otherSlider2Value.html(newStepValue);

    var slider_value = slider.val();
    var selected_gpcr = $('#selected').data('gpcr').replace('"', "").replace('"', '');
    var selected_gprotein = $('#selected').data('gprotein').replace('"', "").replace('"', '');
    var slider1A_value = document.getElementById('slider1A_value').innerHTML;
    var displayName = document.getElementById('displayName').innerHTML;
    var uniq_id = document.getElementById('uniq_id').value.replace('"', "").replace('"', '');
    var path_to_fasta = document.getElementById('path_to_fasta').value.replace('"', "").replace('"', '');

    makeStructure(selected_gpcr, selected_gprotein, slider1A_value, slider_value, uniq_id);
    makeHeatmap(slider1A_value, slider_value, selected_gpcr, selected_gprotein, displayName);
    makeAttentionmap(uniq_id, selected_gpcr, selected_gprotein, displayName)
    setDisplayMenu(path_to_fasta, slider1A_value, slider_value, uniq_id, selected_gpcr, selected_gprotein, displayName);
    makeSequence(selected_gpcr, path_to_fasta, selected_gprotein, slider1A_value, slider_value, uniq_id, displayName);
    makePCA(uniq_id, 'TGF', '33', selected_gpcr, selected_gprotein);

}