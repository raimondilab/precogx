

function addElement (el) {
  Object.assign(el.style, {
    position: "absolute",
    zIndex: 10
  })
  stage.viewer.container.appendChild(el)
}

function createElement (name, properties, style) {
  var el = document.createElement(name)
  el.className = "btn btn-primary btn-sm";
  Object.assign(el, properties)
  Object.assign(el.style, style)
  //Object.assign(el.classList.add, "btn btn-primary btn-sm")

  return el
}

function createElement2 (id, name, properties, style, labelName) {
  var div = document.getElementById("div"+id);
  divID = 'div' + id
  if(typeof(div) == 'undefined' || div == null){
    var div = document.createElement('DIV')
    div.setAttribute("id", divID);
    div.className = "form-check form-switch";
    div.style.font = "bold 15px arial,serif";
    //div.style.top = "120px";
    //div.style.left = "10px";
  }
  Object.assign(div.style, style);


  var el = document.getElementById(id);
  if(typeof(el) == 'undefined' || el == null){
    var el = document.createElement(name)
    el.className = "form-check-input";
  }
  el.setAttribute("id", id);
  Object.assign(el, properties)
  //Object.assign(el.style, style)

  var label = document.getElementById("label"+id);
  if(typeof(label) == 'undefined' || label == null){
    var label = document.createElement("LABEL");
    label.className ="form-check-label"
    //label.style.backgroundColor = "red";
  }
  label.setAttribute("id", "label"+id);
  label.htmlFor = id;
  label.innerHTML = labelName

  div.appendChild(el)
  div.appendChild(label)
  return div
}

// Function to take GPCR/G-protein chains and pdbid and render it to the
// Structure panel
function showStructure(uniq_id, gpcr, chainGPCR, chainGPROT, pdbid, positions, num_contacts, pair_positions) {
  document.getElementById("PDBbutton").innerHTML = pdbid;
  var cutoff = 0.0;
  $.ajax({
    url:"/convertPositionsBW2PDB", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({pdbID: pdbid, positions: positions, pair_positions: pair_positions, num_contacts: num_contacts, gpcr: gpcr, uniq_id: uniq_id}),
    success: function(response){
				//console.log(response);
        mutation_position = response['mutation_position'];
        modified_positions = response['modified_positions'];
        modified_positions_labels = response['modified_positions_labels'];
        modified_num_contacts = response['modified_num_contacts'];
        //alert(modified_num_contacts);
        modified_pair_positions = response['modified_pair_positions'];
        modified_positions_array = modified_positions.split('_');
        modified_positions_labels_array = modified_positions_labels.split('_');
        modified_num_contacts_array = modified_num_contacts.split('_');
        //alert(response['modified_positions_labels']);
        // Define format of selection of given positions (contacts)
        /*
        selection = '-';
        if (modified_positions_array.length) {
          for (var i = 0; i < modified_positions_array.length; i++) {
            //alert (modified_positions[i]);
            selection += '(' + modified_positions_array[i] + " and .CA and :" + chainGPCR + ') or ';
          }
        }
        */

        // Define format of selection of mutation position (if any)
        if (modified_pair_positions == '') {
          modified_pair_positions_array = []
        }
        else {
          modified_pair_positions_array = modified_pair_positions.split('_');
        }
        selectionDistanceEnrichment = [];
        selectionDistanceDepletion = [];
        if (modified_pair_positions_array.length) {
          for (var i = 0; i < modified_pair_positions_array.length; i++) {
            var row = [];
            var pos1 = modified_pair_positions_array[i].split(':')[0];
            var pos2 = modified_pair_positions_array[i].split(':')[1];
            var score = Number(modified_pair_positions_array[i].split(':')[2]);
            if (score >= 0.0) {
              //alert(score+'up');
              row.push(pos1+'.CA and :'+chainGPCR);
              row.push(pos2+'.CA and :'+chainGPCR);
              //alert(row);
              selectionDistanceEnrichment.push(row);
            }
            else {
              //alert(score+'up');
              row.push(pos1+'.CA and :'+chainGPCR);
              row.push(pos2+'.CA and :'+chainGPCR);
              //alert(row);
              selectionDistanceDepletion.push(row);
            }
          }
        }
        stage.removeAllComponents();
        stage.setParameters({backgroundColor: "white"});
        //alert(modified_positions_array);


        stage.loadFile("rcsb://"+pdbid+".cif").then(function (o) {
                  //o.autoView();

                  var pa = o.structure.getPrincipalAxes();
                  stage.animationControls.rotate(pa.getRotationQuaternion(), 1500);

                  var downloadButton = createElement("input", {
                    class:"custom-control-input",
                    type: "button",
                    value: "Download"
                  }, { top: "10px", left: "10px" })
                  downloadButton.onclick = function (e) {
                    stage.makeImage( {
                        factor: 1,
                        antialias: true,
                        trim: false,
                        transparent: false
                    } ).then( function( blob ){
                        NGL.download( blob, "precogx.png" );
                    } );
                  }
                  addElement(downloadButton);

                  var centerAllButton = createElement("input", {
                    type: "button",
                    value: "Center"
                  }, { top: "10px", left: "140px" })
                  centerAllButton.onclick = function (e) {
                    stage.autoView()
                  }
                  addElement(centerAllButton)

                  //alert(mutation_position);

                  var centerMutationButton = createElement("input", {
                    type: "button",
                    value: "Center mutation"
                  }, { top: "10px", left: "250px"})
                  centerMutationButton.onclick = function (e) {
                    o.autoView(String(mutation_position)+':'+chainGPCR+'.CA')
                  }
                  addElement(centerMutationButton)

                  var enrichButton = createElement2("enrich", "input", {
                    type: "checkbox",
                    checked: true,
                    //value: "Enrichment/Depletion"
                  }, { top: "80px", left: "10px" }, "Enriched contact pairs")
                  addElement(enrichButton);
                  document.getElementById('enrich').onclick = function() {
                    //alert(document.getElementById("el1").checked)
                    //DEPLETION.toggleVisibility();
                    ENRICHMENT.toggleVisibility()
                  }

                  var depleteButton = createElement2("deplete", "input", {
                    type: "checkbox",
                    checked: true,
                    //value: "Enrichment/Depletion"
                  }, { top: "110px", left: "10px" }, "Depleted contact pairs")
                  addElement(depleteButton);
                  document.getElementById('deplete').onclick = function() {
                    //alert(document.getElementById("el1").checked)
                    DEPLETION.toggleVisibility();
                    //ENRICHMENT.toggleVisibility()
                  }

                  var labelButton = createElement2("labeling", "input", {
                    type: "checkbox",
                    checked: true,
                    //value: "Enrichment/Depletion"
                  }, { top: "140px", left: "10px" }, "Labels")
                  addElement(labelButton);
                  document.getElementById('labeling').onclick = function() {
                    //alert(document.getElementById("el1").checked)
                    LABELS.toggleVisibility();
                    //ENRICHMENT.toggleVisibility()
                  }

                  var screenButton = createElement2("fullscreen", "input", {
                    type: "checkbox",
                    checked: false
                    //value: "Enrichment/Depletion"
                  }, { top: "170px", left: "10px" }, "Fullscreen")
                  addElement(screenButton);
                  document.getElementById('fullscreen').onclick = function() {
                    //alert(document.getElementById("el1").checked)
                    stage.toggleFullscreen();
                    //ENRICHMENT.toggleVisibility()
                  }


                  var GPCR = o.addRepresentation("cartoon", {
                      sele: ":"+chainGPCR,
                      name: chainGPCR,
                      color: "silver",
                      //color: schemeId,
                  });
                  var GPROTEIN = o.addRepresentation("cartoon", {
                      sele: ":"+chainGPROT,
                      name: chainGPROT,
                      color: "skyblue",
                      //color: schemeId,
                  });


                  var bg_color = { }
                  var sele = new NGL.Selection("")
                  var view = o.structure.getView(sele)
                  var labelText = { }
                  for (var i = 0; i < modified_positions_array.length; i++) {
                      var seleString = modified_positions_array[i];
                      sele.setString(seleString+':'+chainGPCR+'.CA')
                      var atomIndex = view.getAtomIndices()[0];
                      //alert(atomIndex);
                      if (atomIndex !== undefined) {
                          labelText[atomIndex] = modified_positions_labels_array[i];
                      }
                  }

                  var LABELS = o.addRepresentation("label", {
                    labelType: "text",
                    labelText: labelText,
                    fontWeight: 'normal',
                    //sdf: true,
                    color: "black",
                    xOffset: 1.5,
                    fixedSize: 0.1,
                    //showBackground: true,
                    //backgroundColor: bg_color,
                    //borderColor: 'blue',
                  })

                  /*
                  var mutationText = { };
                  var sele2 = new NGL.Selection("");
                  var view2 = o.structure.getView(sele2);
                  var seleString2 = mutation_position;
                  sele2.setString(seleString2+':'+chainGPCR+'.CA');
                  var atomIndex2 = view2.getAtomIndices()[0];
                  alert(atomIndex2);
                  if (atomIndex2 !== undefined) {
                      mutationText[atomIndex2] = response['mutation_position_label'];
                  }

                  alert(mutationText[atomIndex2]);
                  var MUTLABELS = o.addRepresentation("label", {
                    labelType: "text",
                    labelText: mutationText,
                    color: "purple",
                    xOffset: 1,
                    //showBackground: true,
                    //backgroundColor: bg_color,
                    //borderColor: 'blue',
                  });
                  */


                  if (modified_positions_array.length) {
                    for (var i = 0; i < modified_positions_array.length; i++) {
                      o.addRepresentation("ball+stick", {
                          //sele: selection,
                          sele: modified_positions_array[i] + " and .CA and :" + chainGPCR,
                          name: 'extra',
                          //radius: '0.5',
                          radius: modified_num_contacts_array[i],
                          color: "khaki"
                          //color: schemeId,
                      });
                    }
                  }

                  o.addRepresentation("ball+stick", {
                      sele: mutation_position+ " and .CA and :" + chainGPCR,
                      name: 'extra',
                      radius: '0.5',
                      color: "purple"
                      //color: schemeId,
                  });


                  if (modified_pair_positions_array != []) {
                    var ENRICHMENT = o.addRepresentation( "distance", {
                        //atomPair: [ [ "280.CA and :R", "325.CA and :R" ] ],
                        atomPair: selectionDistanceEnrichment,
                        color: "green",
                        labelSize: 0
                      }
                    );
                    var DEPLETION = o.addRepresentation( "distance", {
                        //atomPair: [ [ "280.CA and :R", "325.CA and :R" ] ],
                        atomPair: selectionDistanceDepletion,
                        color: "red",
                        labelSize: 0
                      }
                    );
                  }
                  o.autoView(':'+chainGPCR);
                  //o.autoView("49:R.CA");
                  //o.autoView(mutation_position+ " and .CA and :" + chainGPCR);
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
function resetPDBlist(uniq_id, gpcr, ordered_pdbs, positions, pair_positions, num_contacts) {
  //alert(positions);
  var new_options = '';
  var new_options = "<input class=\"form-control\" type=\"text\" id=\"PDBsearch\" placeholder=\"Search..\">";
  //alert(gpcr);
  for (var i = 0; i < ordered_pdbs.length; i++) {
    var x = ordered_pdbs[i].split('_');
    var pdbid = x[0];
    var chainGPCR = x[1];
    var chainGPROT = x[2];
    //var positions = '1.5';
    new_options += "<li><a class=\"dropdown-item\" onClick=\"showStructure(\'"+uniq_id+"\',\'"+gpcr+"\',\'"+chainGPCR+"\',\'"+chainGPROT+"\',\'"+pdbid+"\',\'"+positions+"\',\'"+num_contacts+"\',\'"+pair_positions+"\')\">"+pdbid+"</a></li>";
    // return only the first values as default to display
    if (i == 0) {
      var first_values = [chainGPCR, chainGPROT, pdbid];
      showStructure(uniq_id, gpcr, first_values[0], first_values[1], first_values[2], positions, num_contacts, pair_positions);
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
        //alert(response['num_contacts']);
        //alert(response['try'][0]+'--'+response['try']);
        resetPDBlist(uniq_id, gpcr, response['ordered_pdbs'], response['positions'], response['pair_positions'], response['num_contacts']);

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
