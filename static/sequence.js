var seq = new Sequence('MALWMRLLPLLALLALWGPGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN');
// Render the sequence with or without rendering options
// (Check the interactive documentation)
seq.render('#sequence-viewer', {
            'showLineNumbers': true,
            'wrapAminoAcids': true,
            'charsPerLine': 70,
            'toolbar': true,
            'title' : "GPCR",
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
