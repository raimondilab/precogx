// Function to read the GPCR name from the Results table and fetch it's FASTA
// and render it in the Sequence Panel
function makeSequence(gpcr) {
  const xhttp = new XMLHttpRequest();
  xhttp.onload = function() {
    var fasta = this.responseText;
    fasta = fasta.split('>');
    //Ignore i = 0 as it's empty
    for (var i = 1; i < fasta.length; i++) {
      var x = fasta[i].split('\n');
      // Check if the first line has the GPCR name; If yes, selec the sequence
      if (x[0] == gpcr) {
        var sequence = '';
        for (var j = 1; j < x.length; j++) {
          sequence += x[j]
        }
      }

    }
    //alert(sequence);
    //var seq = new Sequence('MALWMRLLPLLALLALWGPGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN');
    var seq = new Sequence(sequence);
    seq.render('#sequence-viewer', {
                'showLineNumbers': true,
                'wrapAminoAcids': true,
                'charsPerLine': 60,
                'toolbar': true,
                'title' : gpcr,
                'search': true,
                'header' : {
                    display:true,
                    searchInTitle :false,
                    unit: "Char",
                    showCpl: true,
                    badgeWithUnit : false
                }
            });
    //highlight randomly chosen region in the sequence
    seq.selection(21,33,"green");
    //document.getElementById("pdblist").innerHTML = new_options;
  }
  xhttp.open("GET", "static/OL820/temp/new_fasta_file.txt");
  xhttp.send();
}
