function resetPDBlist(selected_gpcr) {
  /*
  $.ajax({
    url:"/help", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({pdbid: pdbid, chainGPCR: chainGPCR, chainGPROT: chainGPROT, gpcr: gpcr}),
    success: function(response){
				console.log(response);
        //alert(response['status']);
			},
			error: function(error){
				console.log(error);
			}
    });
    */
  var dropDown = document.getElementById('PDBbutton');
  /*
  dropDown.addEventListener("show.bs.dropdown", function () {
    alert("Button is clicked");
    //var selected_gpcr = 'DRD1';
    resetPDBlist(selected_gpcr);
  });
  */
  alert('selected_gpcr '+selected_gpcr);
  //var selected_gpcr = 'hello';
  const xhttp = new XMLHttpRequest();
  xhttp.onload = function() {
    var options = this.responseText;
    options = options.split('\n');
    var new_options = '';
    for (var i = 0; i < options.length; i++) {
      var x = options[i].split(' ');
      var pdbid = x[0];
      var chainGPCR = x[1];
      var chainGPROT = x[2];
      new_options += "<li><a class=\"dropdown-item\" onClick=\"makeStructure(\'"+pdbid+"\',\'"+chainGPCR+"\',\'"+chainGPROT+"\',\'"+selected_gpcr+"\')\">" + pdbid+"</a></li>";
    }
    //alert(options);
    document.getElementById("pdblist").innerHTML = new_options;
  }
  xhttp.open("GET", "static/pdblist.txt");
  xhttp.send();
}
