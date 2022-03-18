// Function to read the GPCR name from the Results table and fetch it's FASTA
// and render it in the Sequence Panel
function makeSequence(gpcr, path_to_fasta, gprotein, cutoff, distance, uniq_id) {
  //alert(gpcr);
  $.ajax({
    url:"/fetchContactsSequence", //the page containing python script
    type: "post", //request type,
    dataType: 'json',
    data: JSON.stringify({gpcr: gpcr, gprotein: gprotein, cutoff: cutoff, distance: distance, path_to_fasta: path_to_fasta, uniq_id: uniq_id}),
    success: function(response){
				console.log(response);
        sequence = response['sequence']
        seq_positions = response['seq_positions']
        bw_positions = response['bw_positions']
        var variant = gpcr.split('_')[1];
        var numberPattern = /\d+/g;
        var variantPosition = Number(variant.match( numberPattern ));
        /*
        if (variantPosition != null && seq_positions.includes(Number(variantPosition)) == false) {
          //alert(variantPosition);
          seq_positions.push(variantPosition);
          seq_positions.sort(function(a, b) {
                                return a - b;
                              });
          }
        seq_positions = [... new Set(seq_positions)];
        */
        //alert(seq_positions);
        //alert(variantPosition);
        var seq = new Sequence(sequence);
        seq.render('#sequence-viewer', {
                    'showLineNumbers': true,
                    'wrapAminoAcids': true,
                    'charsPerLine': 70,
                    'toolbar': true,
                    'title' : gpcr,
                    'search': true,
                    'badge': false,
                    'header' : {
                        display:true,
                        searchInTitle :false,
                        unit: "Char",
                        showCpl: true,
                        badgeWithUnit : false
                    }
                });

        //highlight randomly chosen region in the sequence
        if (seq_positions.length > 0) {
          var example = [];
          for (var i = 0; i < seq_positions.length; i++) {
            if (Number(variantPosition) === seq_positions[i]) {
              example.push({start: Number(seq_positions[i]-1), end: Number(seq_positions[i]), color: "black", underscore: false, bgcolor: "violet", tooltip: 'BW: '+bw_positions[i]});
            }
            else {
              //example.push(seq_positions[j]);
              example.push({start: Number(seq_positions[i]-1), end: Number(seq_positions[i]), color: "black", underscore: false, bgcolor: "khaki", tooltip: 'BW: '+bw_positions[i]});
            }
          }
          //seq.selection(35,43,"blue");
          seq.coverage(example);

          //alert(variant);

          if (variant != 'WT') {
            var exampleLegend = [
                {name: "Mutation: "+variant, color: "violet", underscore: false},
                {name: "Predicted contact positions for "+gprotein+" coupling group", color: "khaki", underscore: false},
                {name: "Hover to view BW annotation"}
                ];
          }
          else {
            var exampleLegend = [
                {name: "Predicted contact positions for "+gprotein+" coupling group", color: "khaki", underscore: false},
                {name: "Hover to view BW annotations"}
                ];
          }
          seq.addLegend(exampleLegend);
        }

			},
			error: function(error){
				console.log(error);
			}
    });

}
