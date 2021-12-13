function makePCA(uniq_id, assay, pca_type, gpcr, gprotein) {
  //showPCA(uniq_id, assay, gpcr, gprotein);
  var gpcr = gpcr;
  var gprotein = gprotein;
  //alert(pca_type+'----'+gpcr+'----'+gprotein);
  $.ajax({
    url:"/fetchPCA", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({uniq_id: uniq_id, assay: assay, pca_type: pca_type, gpcr: gpcr, gprotein: gprotein}),
    success: function(response){
				console.log(response);
        showPCA(uniq_id, assay, pca_type, gpcr, gprotein);
        var train_grey = {
          x:response['x_train_grey'],
          y:response['y_train_grey'],
          mode: 'markers',
          type: 'scatter',
          name: assay + ' grey',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_grey'],
          marker: { size: 12, color: 'lightgrey' }
        };

        var train_coupling = {
          x:response['x_train_coupling'],
          y:response['y_train_coupling'],
          mode: 'markers',
          type: 'scatter',
          name: assay + ' coupling',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_coupling'],
          marker: { size: 12, color: 'orange' }
        };

        var train_uncoupling = {
          x:response['x_train_uncoupling'],
          y:response['y_train_uncoupling'],
          mode: 'markers',
          type: 'scatter',
          name: assay + ' not coupling',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_uncoupling'],
          marker: { size: 12, color: 'blue' }
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
          var data = [ train_grey, train_coupling, train_uncoupling, test, wt ];
        }
        else {
          var data = [ train_grey, train_coupling, train_uncoupling, test ];
        }

        var layout = {
          xaxis: {
            range: [ response['minX'], response['maxX'] ]
          },
          yaxis: {
            range: [ response['minY'], response['maxY'] ]
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

function showPCA(uniq_id, assay, pca_type, gpcr, gprotein) {
  //alert(pca_type);
  //alert(gpcr+'--'+gprotein);
  var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'Shedding\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Shedding</a></li>";
  options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'ebBRET\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">ebBRET</a></li>";
  options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'IUPHAR\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">IUPHAR</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
  //alert(options);
  document.getElementById("AssayList").innerHTML = options;
  document.getElementById("AssayButton").innerHTML = assay;
  //showPCA(uniq_id, assay, gpcr, gprotein);

  var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
  options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
  document.getElementById("PCAList").innerHTML = options;
  document.getElementById("PCAButton").innerHTML = pca_type;
}
