function makeAttentionmap(uniq_id, gpcr, gprotein, chosen) {
  var gpcr = gpcr;
  var gprotein = gprotein;
  $.ajax({
    url:"/fetchAttentionMap", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({gpcr: gpcr, gprotein: gprotein, uniq_id: uniq_id}),
    success: function(response){
				console.log(response);
        togglePanel(chosen);
        var data = [
          {
            //z: [[1, 20, 30], [20, 1, 60], [30, 60, 1]],
            z: response['zaxis'],
            //zmin: Number(response['fetch_contactsMin']),
            //zmax: Number(response['fetch_contactsMax']),
            //zmin: -1.0,
            //zmax: 1.0,
            x: response['xaxis'],
            y: response['yaxis'],
            type: 'heatmap',
            colorscale: 'Reds',
            /*
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
            */
            hoverongaps: false
          }
        ];

        var layout = {
              autosize: true,
              title: 'Attention Map of <b>'+gpcr+'</b> with the <b>'+gprotein+'</b> coupling group',
              font: {
                size: 10
               }

        };

        var icon1 = {
            'width': 500,
            'height': 600

}
        var config = {responsive: true,
        displaylogo: false

        }


        Plotly.newPlot('attentionMap', data, layout, config);

			},
			error: function(error){
				console.log(error);
			}
    });
}

function togglePanel(chosen) {
  // get the clock
  var contactPanel = document.getElementById('myDiv');
  var sequencePanel = document.getElementById('sequence-viewer');
  var attentionPanel = document.getElementById('attentionMap');
  var panelHeading = document.getElementById('panelHeading');
  var slider1Div = document.getElementById('slider1Div');

  // get the current value of the clock's display property
  var contactDisplaySetting = contactPanel.style.display;
  var sequenceDisplaySetting = sequencePanel.style.display;
  var attentionDisplaySetting = attentionMap.style.display;
  var slider1DisplaySetting = slider1Div.style.display;

  // also get the clock button, so we can change what it says
  //var toggleButton = document.getElementById('toggleButton');

  // now toggle the clock and the button text, depending on current state
  if (chosen == 'Contact Map') {
    // clock is visible. hide it
    contactPanel.style.display = 'block';
    sequencePanel.style.display = 'none';
    attentionPanel.style.display = 'none';
    slider1Div.style.display = 'block';
    //makeHeatmap();
    // change button text
    //toggleButton.innerHTML = 'Map the contacts on sequence';
    //panelHeading.innerHTML = 'Differential Contact Map';
    displayName.innerHTML = chosen;
  }
  else if (chosen == 'Attention Map') {
    // clock is hidden. show it
    contactPanel.style.display = 'none';
    sequencePanel.style.display = 'none';
    attentionPanel.style.display = 'block';
    slider1Div.style.display = 'none';
    // change button text
    //toggleButton.innerHTML = 'Show the contact pairs';
    //panelHeading.innerHTML = 'Attention Map';
    displayName.innerHTML = chosen;
  }
  else {
    // clock is hidden. show it
    contactPanel.style.display = 'none';
    sequencePanel.style.display = 'block';
    attentionPanel.style.display = 'none';
    slider1Div.style.display = 'block';
    // change button text
    //toggleButton.innerHTML = 'Show the contact pairs';
    //panelHeading.innerHTML = 'Sequence';
    displayName.innerHTML = chosen;
  }
}

function setDisplayMenu(path_to_fasta, cutoff, distance, uniq_id, gpcr, gprotein) {
  var new_options = '';
  new_options += "<li><a class=\"dropdown-item\" onClick=\"makeHeatmap(\'"+cutoff+"\',\'"+distance+"\',\'"+gpcr+"\',\'"+gprotein+"\',\'Contact Map\')\">Contact Map</a></li>";
  new_options += "<li><a class=\"dropdown-item\" onClick=\"makeAttentionmap(\'"+uniq_id+"\',\'"+gpcr+"\',\'"+gprotein+"\',\'Attention Map\')\">Attention Map</a></li>";
  new_options += "<li><a class=\"dropdown-item\" onClick=\"makeSequence(\'"+gpcr+"\',\'"+path_to_fasta+"\',\'"+gprotein+"\',\'"+cutoff+"\',\'"+distance+"\',\'"+uniq_id+"\',\'Sequence\')\">Sequence</a></li>";

  //alert(document.getElementById("pdblist").innerHTML);
  //document.getElementById("displayList").innerHTML = new_options;
  document.getElementById("displayList").innerHTML = new_options;
  }
