function makeStructure(pdbid, chainGPCR, chainGPROT) {
  //var id = pdbid;
  stage.removeAllComponents();
  stage.setParameters( { backgroundColor: "white"} );
  //stage.loadFile("static/"+id+".cif", {defaultRepresentation: true});
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

function makeStructure2(val) {
  //alert(id);
  //var stage;
  var id = val;
  //var pdb_id = id;
  //stage = new NGL.Stage("viewport");
  //alert(id);
  stage.removeAllComponents();
  stage.setParameters( { backgroundColor: "white"} );
  stage.loadFile("static/"+id+".cif").then(function (o) {
            o.autoView();
            o.addRepresentation("cartoon", {
                sele: ":R",
                name: "R",

                color: "blue",
                //color: schemeId,
            });
            o.addRepresentation("cartoon", {
                sele: ":A",
                name: "A",

                color: "green",
                //color: schemeId,
            });
          });
}
