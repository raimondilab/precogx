function makePCA2(uniq_id, assay, pca_type, gpcr, gprotein) {
  //showPCA(uniq_id, assay, gpcr, gprotein);
  var gpcr = gpcr;
  var gprotein = gprotein;
  //alert(assay+'----'+gpcr+'----'+gprotein);
  //var assay = showPCA(uniq_id, assay, pca_type, gpcr, gprotein);
  $.ajax({
    url:"/fetchPCA2", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({uniq_id: uniq_id, assay: assay, pca_type: pca_type, gpcr: gpcr, gprotein: gprotein}),
    success: function(response){
				console.log(response);
        //alert(assay);
        //alert(response['y_classA']);
        var assay = response['assay'];
        showPCA(uniq_id, assay, pca_type, gpcr, gprotein);
        var other = {
          x:response['x_other'],
          y:response['y_other'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Other',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['other'],
          marker: { size: 12, color: 'lightgrey' }
        };

        var classA = {
          x:response['x_classA'],
          y:response['y_classA'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Class A',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['classA'],
          marker: { size: 12, color: 'lightgreen' }
        };

        var classB1 = {
          x:response['x_classB1'],
          y:response['y_classB1'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Class B1',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['classB1'],
          marker: { size: 12, color: 'steelblue' }
        };

        var classB2 = {
          x:response['x_classB2'],
          y:response['y_classB2'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Class B2',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['classB2'],
          marker: { size: 12, color: 'cyan' }
        };

        var classC = {
          x:response['x_classC'],
          y:response['y_classC'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Class C',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['classC'],
          marker: { size: 12, color: 'forestgreen' }
        };

        var taste = {
          x:response['x_taste'],
          y:response['y_taste'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Taste',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['taste'],
          marker: { size: 12, color: 'darkgrey' }
        };

        var frizzeled = {
          x:response['x_frizzeled'],
          y:response['y_frizzeled'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Frizzeled',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['frizzeled'],
          marker: { size: 12, color: 'khaki' }
        };

        if (response['x_wt'] != '-') {
          //var testColor = "fuchsia";
          //var testColor = "#0000FF";
          var testColor = "#FF0000";
        }
        else {
          //var testColor = "darkorange";
          //var testColor = "#FF0000";
          var testColor = "#0000FF";
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
            //marker: { size: 12, color: "darkorange" }
            marker: { size: 12, color: "#0000FF" }
          };
          var data = [ other, classA, classB1, classB2, classC, frizzeled, taste, test, wt ];
          //var data = [ classA, test, wt ];
        }
        else {
          var data = [ other, classA, classB1, classB2, classC, frizzeled, taste, test ];
          //var data = [ classA, test ];
        }

        var layout = {
          xaxis: {
            range: [ response['minX'], response['maxX'] ]
          },
          yaxis: {
            range: [ response['minY'], response['maxY'] ]
          },
          title: gprotein+'--'
        };

        var config = {responsive: true,
        displaylogo: false}

        Plotly.newPlot('myDiv2', data, layout, config);
			},
			error: function(error){
				console.log(error);
			}
    });
}

function makePCA(uniq_id, assay, pca_type, gpcr, gprotein) {
  //showPCA(uniq_id, assay, gpcr, gprotein);
  var gpcr = gpcr;
  var gprotein = gprotein;
  //alert(assay+'----'+gpcr+'----'+gprotein);
  //var assay = showPCA(uniq_id, assay, pca_type, gpcr, gprotein);
  $.ajax({
    url:"/fetchPCA", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({uniq_id: uniq_id, assay: assay, pca_type: pca_type, gpcr: gpcr, gprotein: gprotein}),
    success: function(response){
				console.log(response);
        //alert(response['score_coupling']);
        var assay = response['assay'];
        //alert(assay);
        showPCA(uniq_id, assay, pca_type, gpcr, gprotein);
        var train_grey = {
          x:response['x_train_grey'],
          y:response['y_train_grey'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' grey',
          name: 'Not annotated',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_grey'],
          marker: { size: 12, color: 'lightgrey' }
        };

        var train_uncoupling = {
          x:response['x_train_uncoupling'],
          y:response['y_train_uncoupling'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' not coupling',
          name: 'Not coupled',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_uncoupling'],
          marker: { size: 12, color: response['score_uncoupling'] }
        };

        var train_coupling = {
          x:response['x_train_coupling'],
          y:response['y_train_coupling'],
          mode: 'markers',
          type: 'scatter',
          //name: assay + ' coupling',
          name: 'Coupled',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_coupling'],
          marker: { size: 12, color: response['score_coupling'] }
        };

        if (response['x_wt'] != '-') {
          //var testColor = "fuchsia";
          var testColor = "#FF0000";
        }
        else {
          //var testColor = "darkorange";
          var testColor = "#0000FF";
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
            marker: { size: 12, color: "#0000FF" }
          };
          var data = [ train_grey, train_uncoupling, train_coupling, test, wt ];
        }
        else {
          var data = [ train_grey, train_uncoupling, train_coupling, test ];
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

        var config = {responsive: true}

        Plotly.newPlot('myDiv2', data, layout, config);
			},
			error: function(error){
				console.log(error);
			}
    });
}

function showPCA(uniq_id, assay, pca_type, gpcr, gprotein) {
  //alert(gpcr+'--'+gprotein);
  var gemta = ['GNAS', 'GNAI1', 'GNAI2', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2'];
  var tgf = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15'];
  var both = ['GNAS', 'GNAI1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GoA'];
  var options = '';

  if (gprotein.includes("Barr")) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GEMTA\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GEMTA</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'STRING\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">STRING</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'Class\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Class</a></li>";
  }
  else if (both.includes(gprotein)) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'TGF\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">TGF</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GEMTA\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GEMTA</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GtoPdb\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GtoPdb</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'Class\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Class</a></li>";
  }
  else if (tgf.includes(gprotein)) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'TGF\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">TGF</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GtoPdb\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GtoPdb</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'Class\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Class</a></li>";
  }
  else {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GEMTA\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GEMTA</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'GtoPdb\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">GtoPdb</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'Class\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Class</a></li>";
  }

  document.getElementById("AssayButton").innerHTML = assay;
  document.getElementById("AssayList").innerHTML = options;

  /*
  if (assay.includes("Class")) {
    var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
    for (let i = 0; i < 34; i++) {
      options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'"+i+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">"+i+"</a></li>";
    }
  }
  else {
    var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
    for (let i = 0; i < 34; i++) {
      options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'"+i+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">"+i+"</a></li>";
    }
  }
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
  document.getElementById("PCAList").innerHTML = options;
  document.getElementById("PCAButton").innerHTML = pca_type;
  */

  $.ajax({
    url:"/assignOptions", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({uniq_id: uniq_id, assay: assay, pca_type: pca_type, gpcr: gpcr, gprotein: gprotein}),
    success: function(response){
				console.log(response);
        //alert(assay);
        //alert(response['score_coupling']);
        var layers = response['layers'];
        if (assay.includes("Class")) {
          var options = '';
          //var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
          //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
          for (let i = 0; i < layers.length; i++) {
            options += "<li><a class=\"dropdown-item\" onClick=\"makePCA2(\'"+uniq_id+"\',\'"+assay+"\'"+",\'"+i+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">"+layers[i]+"</a></li>";
          }
        }
        else {
          var options = '';
          //var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
          //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
          for (let i = 0; i < layers.length; i++) {
            options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'"+i+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">"+layers[i]+"</a></li>";
          }
        }
        //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
        //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
        document.getElementById("PCAList").innerHTML = options;
        document.getElementById("PCAButton").innerHTML = pca_type;

			},
			error: function(error){
				console.log(error);
			}
    });
}
