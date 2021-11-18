// Function to take GPCR/G-protein chains and pdbid and render it to the
// Structure panel
function showStructure(chainGPCR, chainGPROT, pdbid) {
  document.getElementById("PDBbutton").innerHTML = pdbid;
  //chainGPCR = 'R';
  //chainGPROT = 'A';
  //pdbid = '3sn6';
  //var id = pdbid;
  //alert(pdbid+' structure');
  stage.removeAllComponents();
  stage.setParameters({backgroundColor: "white"});
  //stage.loadFile("static/"+pdbid+".cif", {defaultRepresentation: true});
  //stage.loadFile("static/"+pdbid+".cif").then(function (o) {
  stage.loadFile("rcsb://"+pdbid+".cif").then(function (o) {
            //o.autoView();
            o.addRepresentation("cartoon", {
                sele: ":"+chainGPCR,
                name: chainGPCR,
                color: "blue",
                //color: schemeId,
            });
            o.addRepresentation("cartoon", {
                sele: ":"+chainGPROT,
                name: chainGPROT,
                color: "green",
                //color: schemeId,
            });
            o.addRepresentation("ball+stick", {
                sele: "1-100 and .CA and :"+chainGPCR,
                name: 'extra',
                radius: '1',
                color: "red",
                //color: schemeId,
            });
            o.autoView();
            //o.center({sele:'1-50',});
          });
}

// Function to take GPCR and orde PDB list as input,
// and insert in the pdblist ID
function resetPDBlist(gpcr, ordered_pdbs) {
  var new_options = '';
  var new_options = "<input class=\"form-control\" type=\"text\" id=\"PDBsearch\" placeholder=\"Search..\">";
  for (var i = 0; i < ordered_pdbs.length; i++) {
    var x = ordered_pdbs[i].split('_');
    var pdbid = x[0];
    var chainGPCR = x[1];
    var chainGPROT = x[2];
    new_options += "<li><a class=\"dropdown-item\" onClick=\"showStructure(\'"+chainGPCR+"\',\'"+chainGPROT+"\',\'"+pdbid+"\')\">" + pdbid+"</a></li>";
    // return only the first values as default to display
    if (i == 0) {
      var first_values = [chainGPCR, chainGPROT, pdbid];
      showStructure(first_values[0], first_values[1], first_values[2]);
    }
  }
  //alert(new_options);
  document.getElementById("pdblist").innerHTML = new_options;
}

// Function to take GPCR as input and make new PDB ordered list
function makeStructure(gpcr) {
  $.ajax({
    url:"/order_pdbs", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    //data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    data: JSON.stringify({gpcr: gpcr}),
    success: function(response){
				console.log(response);
        //alert(response['ordered_pdbs']);
        resetPDBlist(gpcr, response['ordered_pdbs']);
        //showStructure(first_values[0], first_values[1], first_values[2]);
        //alert(first_values);

        //set search filter
        $("#PDBsearch").on("keyup", function() {
          var value = $(this).val().toLowerCase();
          //alert(value);
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
