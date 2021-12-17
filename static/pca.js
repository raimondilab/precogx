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
        //alert(assay);
<<<<<<< HEAD
=======
        //alert(response['score_coupling']);
>>>>>>> refs/remotes/origin/main
        var assay = response['assay'];
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
          marker: { size: 12, color: response['score_coupling'] }
        };

        var train_uncoupling = {
          x:response['x_train_uncoupling'],
          y:response['y_train_uncoupling'],
          mode: 'markers',
          type: 'scatter',
          name: assay + ' not coupling',
          //text: ['A-1', 'A-2', 'A-3', 'A-4', 'A-5'],
          text: response['genes_to_consider_uncoupling'],
          marker: { size: 12, color: response['score_uncoupling'] }
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
            marker: { size: 12, color: "darkblue" }
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
  //alert(gpcr+'--'+gprotein);
  var ebBRET = ['GNAS', 'GNAI1', 'GNAI2', 'GoA', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2'];
  var shedding = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAO1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15'];
  var both = ['GNAS', 'GNAI1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15'];
  var options = '';

  if (gprotein.includes("Barr")) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'ebBRET\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">ebBRET</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'STRING\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">STRING</a></li>";
  }
  else if (both.includes(gprotein)) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'Shedding\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Shedding</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'ebBRET\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">ebBRET</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'IUPHAR\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">IUPHAR</a></li>";
  }
  else if (shedding.includes(gprotein)) {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'Shedding\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">Shedding</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'IUPHAR\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">IUPHAR</a></li>";
  }
  else {
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'ebBRET\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">ebBRET</a></li>";
    options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'IUPHAR\'"+",\'"+pca_type+"\',\'"+gpcr+"\',\'"+gprotein+"\')\">IUPHAR</a></li>";
  }

  document.getElementById("AssayButton").innerHTML = assay;
  document.getElementById("AssayList").innerHTML = options;

  var options = "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'GPCRome\',\'"+gpcr+"\',\'"+gprotein+"\')\">GPCRome</a></li>";
  options += "<li><a class=\"dropdown-item\" onClick=\"makePCA(\'"+uniq_id+"\',\'"+assay+"\'"+",\'Best PCA\',\'"+gpcr+"\',\'"+gprotein+"\')\">Best PCA</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">Shedding</a></li>";
  //options += "<li><a class=\"dropdown-item\" onClick=\"makePCA()\">ebBRET</a></li>";
  document.getElementById("PCAList").innerHTML = options;
  document.getElementById("PCAButton").innerHTML = pca_type;
}
