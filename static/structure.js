// Function to read the PDB-ID, GPCR chain and G-protein chain from pdblist.txt
// and render it to the Structure Panel
function makeStructure(pdbid, chainGPCR, chainGPROT, gpcr) {
  $.ajax({
    url:"/help", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    success: function(response){
				console.log(response);
        alert(response['status']);
			},
			error: function(error){
				console.log(error);
			}
    });
  //var id = pdbid;
  //alert(gpcr+'structure');
  stage.removeAllComponents();
  stage.setParameters( { backgroundColor: "white"} );
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
