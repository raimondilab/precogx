// Function to take GPCR/G-protein chains and pdbid and render it to the
// Structure panel
function showStructure(uniq_id, gpcr, chainGPCR, chainGPROT, pdbid, positions, pair_positions) {
  //alert (positions);
  document.getElementById("PDBbutton").innerHTML = pdbid;
  var cutoff = 0.0;
  $.ajax({
    url:"/convertPositionsBW2PDB", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({pdbID: pdbid, positions: positions, pair_positions: pair_positions, gpcr: gpcr, uniq_id: uniq_id}),
    success: function(response){
				//console.log(response);
        mutation_position = response['mutation_position'];
        modified_positions = response['modified_positions'];
        modified_pair_positions = response['modified_pair_positions'];
        modified_positions_array = modified_positions.split('_');
        //alert(modified_pair_positions);
        // Define format of selection of given positions (contacts)
        selection = '-';
        if (modified_positions_array.length) {
          for (var i = 0; i < modified_positions_array.length; i++) {
            //alert (modified_positions[i]);
            selection += '(' + modified_positions_array[i] + " and .CA and :" + chainGPCR + ') or ';
          }
        }

        // Define format of selection of mutation position (if any)
        modified_pair_positions_array = modified_pair_positions.split('_');
        selectionDistance = [];
        if (modified_pair_positions_array.length) {
          for (var i = 0; i < modified_pair_positions_array.length; i++) {
            var row = [];
            var pos1 = modified_pair_positions_array[i].split('-')[0];
            var pos2 = modified_pair_positions_array[i].split('-')[1];
            row.push(pos1+'.CA and :'+chainGPCR);
            row.push(pos2+'.CA and :'+chainGPCR);
            //alert(row);
            selectionDistance.push(row);
          }
        }
        stage.removeAllComponents();
        stage.setParameters({backgroundColor: "white"});

        stage.loadFile("rcsb://"+pdbid+".cif").then(function (o) {
                  //o.autoView();
                  o.addRepresentation("cartoon", {
                      sele: ":"+chainGPCR,
                      name: chainGPCR,
                      color: "silver",
                      //color: schemeId,
                  });
                  o.addRepresentation("cartoon", {
                      sele: ":"+chainGPROT,
                      name: chainGPROT,
                      color: "lightgreen",
                      //color: schemeId,
                  });
                  o.addRepresentation("ball+stick", {
                      sele: selection,
                      name: 'extra',
                      radius: '1',
                      color: "khaki",
                      //color: schemeId,
                  });
                  o.addRepresentation("ball+stick", {
                      sele: mutation_position+ " and .CA and :" + chainGPCR,
                      name: 'extra',
                      radius: '1',
                      color: "red",
                      //color: schemeId,
                  });
                  o.addRepresentation( "distance", {
                      //atomPair: [ [ "280.CA and :R", "325.CA and :R" ] ],
                      atomPair: selectionDistance,
                      color: "skyblue"
                  } );
                  o.autoView(':'+chainGPCR);
                  //o.removeAllRepresentations();
                  //o.center({sele:'1-50',});
                });
        //alert(response['modified_positions']);
			},
			error: function(error){
				console.log(error);
        alert('structure error -- probably bcoz contacts unavailable for the given gprotein');
			}

    });
  }

// Function to take GPCR and orde PDB list as input,
// and insert in the pdblist ID
function resetPDBlist(uniq_id, gpcr, ordered_pdbs, positions, pair_positions) {
  var new_options = '';
  var new_options = "<input class=\"form-control\" type=\"text\" id=\"PDBsearch\" placeholder=\"Search..\">";
  //alert(gpcr);
  for (var i = 0; i < ordered_pdbs.length; i++) {
    var x = ordered_pdbs[i].split('_');
    var pdbid = x[0];
    var chainGPCR = x[1];
    var chainGPROT = x[2];
    //var positions = '1.5';
    new_options += "<li><a class=\"dropdown-item\" onClick=\"showStructure(\'"+uniq_id+"\',\'"+gpcr+"\',\'"+chainGPCR+"\',\'"+chainGPROT+"\',\'"+pdbid+"\',\'"+positions+"\',\'"+pair_positions+"\')\">"+pdbid+"</a></li>";
    // return only the first values as default to display
    if (i == 0) {
      var first_values = [chainGPCR, chainGPROT, pdbid];
      showStructure(uniq_id, gpcr, first_values[0], first_values[1], first_values[2], positions, pair_positions);
    }
    //alert('hello');
  }
  document.getElementById("pdblist").innerHTML = new_options;
}

// Function to take GPCR as input and make new PDB ordered list
function makeStructure(gpcr, gprotein, cutoff, uniq_id) {
  //alert(gprotein);
  //var cutoff =0.0;
  $.ajax({
    url:"/fetchContactsPDBStructure", //to fetch contacts and ordered PDB list based on G-protein
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({gpcr: gpcr, gprotein: gprotein, cutoff: cutoff, uniq_id: uniq_id}),
    success: function(response){
				console.log(response);
        //alert(response['ordered_pdbs']);
        resetPDBlist(uniq_id, gpcr, response['ordered_pdbs'], response['positions'], response['pair_positions']);

        //set search filter
        $("#PDBsearch").on("keyup", function() {
          var value = $(this).val().toLowerCase();
          $(".dropdown-menu li").filter(function() {
            $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
          });

        });
			},
			error: function(error){
				console.log(error);
			}
    });

  }
