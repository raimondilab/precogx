function readOutput() {
  //alert('HELOO');
  var dataSet = [];
  const xhttp = new XMLHttpRequest();
  xhttp.onload = function() {
    var options = this.responseText;
    //alert(options);
    //var options = "DRD1\nDRD2\nDRD3\nDRD4";
    options = options.split('\n');
    //var new_options = '<thead><th>GENE</th></thead><tbody>';
    var new_options = '';
    for (var i = 0; i < options.length; i++) {
      if (i == 0) {
        var x = options[i].split('\t');
        new_options += '<thead>\n<tr>\n';
        for (var j = 0; j < x.length; j++) {
          new_options += "<th>"+x[j].toString()+"</th>\n";
        }
      }
      else if (i < options.length-1) {
        var x = options[i].split('\t');
        var row = [];
        for (var j = 0; j < x.length; j++) {
          row.push(x[j].toString())
        }
        dataSet.push(row);
      }
    }
    //alert(dataSet);
    //return(dataSet);
    dataT.ajax.reload()
    //document.getElementById("example").innerHTML = new_options;
  }
  xhttp.open("GET", "static/OL820/out.tsv");
  xhttp.send();
  //alert (dataSet);
  var myJsonString = JSON.stringify(dataSet);
  //alert (myJsonString);
  return (dataSet);
}
