function makeSequence(gpcr) {
  const xhttp = new XMLHttpRequest();
  xhttp.onload = function() {
    var options = this.responseText;
    options = options.split('>');
    var new_options = '';
    //Ignore i = 0 as it's empty
    for (var i = 1; i < options.length; i++) {
      var x = options[i].split('\n');
      if (x[0] == gpcr) {
        //alert ('found');
        var sequence = '';
        for (var j = 1; j < x.length; j++) {
          sequence += x[j]
        }
      }
      //new_options += "<li><a class=\"dropdown-item\" onClick=\"makeStructure(\'"+pdbid+"\',\'"+chainGPCR+"\',\'"+chainGPROT+"\')\">" + pdbid+"</a></li>";
    }
    //alert(sequence);
    var seq = new Sequence('MALWMRLLPLLALLALWGPGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN');
    var seq = new Sequence(sequence);
    seq.render('#sequence-viewer', {
                'showLineNumbers': true,
                'wrapAminoAcids': true,
                'charsPerLine': 70,
                'toolbar': true,
                'title' : gpcr,
                //'search': true,
                'header' : {
                    display:true,
                    searchInTitle :true,
                    unit: "Char",
                    showCpl: true,
                    badgeWithUnit : false
                }
            });
    seq.selection(21,33,"green");
    //fasta(sequence);
    //document.getElementById("pdblist").innerHTML = new_options;
  }
  xhttp.open("GET", "static/OL820/temp/new_fasta_file.txt");
  xhttp.send();
  /*
}
function fasta(sequence) {
  */

}
