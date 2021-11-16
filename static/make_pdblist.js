function loadDoc() {
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
      new_options += "<li><a class=\"dropdown-item\" onClick=\"makeStructure(\'"+pdbid+"\',\'"+chainGPCR+"\',\'"+chainGPROT+"\')\">" + pdbid+"</a></li>";
    }
    //alert(options);
    document.getElementById("pdblist").innerHTML = new_options;
  }
  xhttp.open("GET", "static/help.txt");
  xhttp.send();
}
