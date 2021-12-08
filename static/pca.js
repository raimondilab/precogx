function makePCA(uniq_id, label, gpcr, gprotein) {
  //showPCA(uniq_id, label, gpcr, gprotein);
  var gpcr = gpcr;
  var gprotein = gprotein;
  //alert(gpcr);
  $.ajax({
    url:"/fetchPCA", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({uniq_id: uniq_id, label: label, gpcr: gpcr, gprotein: gprotein}),
    success: function(response){
				console.log(response);
        showPCA(uniq_id, label, gpcr, gprotein);
        var train_coupling = {
          x:response['x_train_coupling'],
          y:response['y_train_coupling'],
          mode: 'markers',
          type: 'scatter',
          name: label + ' coupling',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_coupling'],
          marker: { size: 12 }
        };

        var train_uncoupling = {
          x:response['x_train_uncoupling'],
          y:response['y_train_uncoupling'],
          mode: 'markers',
          type: 'scatter',
          name: label + ' not coupling',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_uncoupling'],
          marker: { size: 12 }
        };

        if (response['x_wt'] != '-') {
          var testColor = "purple";
        }
        else {
          var testColor = "green";
        }
        var test = {
          x: [response['x_test']],
          y: [response['y_test']],
          mode: 'markers',
          type: 'scatter',
          name: gpcr,
          //text: ['B-a', 'B-b', 'B-c', 'B-d', 'B-e'],
          marker: { size: 12, color: testColor }
        };

        if (response['x_wt'] != '-') {
          var wt = {
            x: [response['x_wt']],
            y: [response['y_wt']],
            mode: 'markers',
            type: 'scatter',
            name: gpcr.split('_')[0]+'_WT',
            //text: ['B-a', 'B-b', 'B-c', 'B-d', 'B-e'],
            marker: { size: 12, color: "green" }
          };
          var data = [ train_coupling, train_uncoupling, test, wt ];
        }
        else {
          var data = [ train_coupling, train_uncoupling, test ];
        }

        var layout = {
          xaxis: {
            range: [ response['min'], response['max'] ]
          },
          yaxis: {
            range: [ response['min'], response['max'] ]
          },
          title: gprotein
        };

        Plotly.newPlot('myDiv2', data, layout);
			},
			error: function(error){
				console.log(error);
			}
    });
}

function showPCA(uniq_id, label, gpcr, gprotein) {
  var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'Shedding\',\'"+gpcr+"\',\'"+gprotein+"\')\">Shedding</a></li>";
  options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'ebBRET\',\'"+gpcr+"\',\'"+gprotein+"\')\">ebBRET</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
  document.getElementById("labelList").innerHTML = options;
  //showPCA(uniq_id, label, gpcr, gprotein);
}
