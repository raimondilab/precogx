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
            //zmin: Number(response['fetch_contactsMin']),
            //zmax: Number(response['fetch_contactsMax']),
            zmin: -1.0,
            zmax: 1.0,
            x: response['positions'],
            y: response['positions'],
            type: 'heatmap',
            colorscale: [
              ['0.0', 'rgb(165,0,0)'],
              ['0.111111111111', 'rgb(215,48,0)'],
              ['0.222222222222', 'rgb(244,109,0)'],
              ['0.333333333333', 'rgb(253,174,0)'],
              ['0.444444444444', 'rgb(254,224,0)'],
              ['0.555555555556', 'rgb(224,243,0)'],
              ['0.666666666667', 'rgb(171,217,0)'],
              ['0.777777777778', 'rgb(116,173,0)'],
              ['0.888888888889', 'rgb(69,117,0)'],
              ['1.0', 'rgb(49,54,0)']
            ],
            hoverongaps: false
          }
        ];

        var layout = {
              autosize: false,
              width: 420,
              height: 420
        };

        Plotly.newPlot('myDiv', data, layout, {displayModeBar: true, scrollZoom: true});

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
