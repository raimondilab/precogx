function makeHeatmap(cutoff, gpcr, gprotein) {
  var gpcr = gpcr;
  var gprotein = gprotein;
  $.ajax({
    url:"/fetchContactsHeatmap", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({gpcr: gpcr, gprotein: gprotein, cutoff: cutoff}),
    success: function(response){
				console.log(response);
        var data = [
          {
            //z: [[1, 20, 30], [20, 1, 60], [30, 60, 1]],
            z: response['fetch_contacts'],
            x: response['positions'],
            y: response['positions'],
            type: 'heatmap'
            //hoverongaps: false
          }
        ];
        Plotly.newPlot('myDiv', data);
			},
			error: function(error){
				console.log(error);
			}
    });
}

function togglePanel() {
  // get the clock
  var contactPanel = document.getElementById('myDiv');
  var sequencePanel = document.getElementById('sequence-viewer');
  var panelHeading = document.getElementById('panelHeading');

  // get the current value of the clock's display property
  var contactDisplaySetting = contactPanel.style.display;
  var sequenceDisplaySetting = sequencePanel.style.display;

  // also get the clock button, so we can change what it says
  var toggleButton = document.getElementById('toggleButton');

  // now toggle the clock and the button text, depending on current state
  if (contactDisplaySetting == 'none') {
    // clock is visible. hide it
    contactPanel.style.display = 'block';
    sequencePanel.style.display = 'none';
    //makeHeatmap();
    // change button text
    toggleButton.innerHTML = 'Show sequence';
    panelHeading.innerHTML = 'Contacts';
  }
  else {
    // clock is hidden. show it
    contactPanel.style.display = 'none';
    sequencePanel.style.display = 'block';
    // change button text
    toggleButton.innerHTML = 'Show contacts';
    panelHeading.innerHTML = 'Sequence';
  }
}
