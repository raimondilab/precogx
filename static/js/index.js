$( document ).ready(function() {

 $(document).on('click', '#ex1', function() {

      var seq1 = ">DRD1\n";
      seq1 += "MRTLNTSAMDGTGLVVERDFSVRILTACFLSLLILSTLLGNTLVCAAVIRFRHLRSKVTN\n";
      seq1 += "FFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFGSFCNIWVAFDIMCSTASILNLCVISVD\n";
      seq1 += "RYWAISSPFRYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAKPTSPSDGNATSLA\n";
      seq1 += "ETIDNCDSSLSRTYAISSSVISFYIPVAIMIVTYTRIYRIAQKQIRRIAALERAAVHAKN\n";
      seq1 += "CQTTTGNGKPVECSQPESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCILPFCGS\n";
      seq1 += "GETQPFCIDSNTFDVFVWFGWANSSLNPIIYAFNADFRKAFSTLLGCYRLCPATNNAIET\n";
      seq1 += "VSINNNGAAMFSSHHEPRGSISKECNLVYLIPHAVGSSEDLKKEEAAGIARPLEKLSPAL\n";
      seq1 += "SVILDYDTDVSLEKIQPITQNGQHPT\n";
      seq1 += ">DRD1/N5P\n"
      var seq2 = ">FFAR2\n";
      seq2 += "MLPDWKSSLILMAYIIIFLTGLPANLLALRAFVGRIRQPQPAPVHILLLSLTLADLLLLL\n";
      seq2 += "LLPFKIIEAASNFRWYLPKVVCALTSFGFYSSIYCSTWLLAGISIERYLGVAFPVQYKLS\n";
      seq2 += "RRPLYGVIAALVAWVMSFGHCTIVIIVQYLNTTEQVRSGNEITCYENFTDNQLDVVLPVR\n";
      seq2 += "LELCLVLFFIPMAVTIFCYWRFVWIMLSQPLVGAQRRRRAVGLAVVTLLNFLVCFGPYNV\n";
      seq2 += "SHLVGYHQRKSPWWRSIAVVFSSLNASLDPLLFYFSSSVVRRAFGRGLQVLRNQGSSLLG\n";
      seq2 += "RRGKDTAEGTNEDRGVGQGEGMPSSDFTTE";

      $('#validationTextarea').html(seq1);
    });

    $(document).on('click', '#ex2', function() {
      var seq2 = "DRD1/F61P\n";
      seq2 += "CXCR3/D278A";

      $('#validationTextarea').html(seq2);
    });

      //$('#validationTextarea').html("DRD1/F61P\nCXCR3/D278A");
      $("#running").hide();

      $("#submit").click(function(){
        $("#home").slideUp();

        var element = document.getElementById("footer");
        element.classList.add("footer");

        $("#running").show(1000);
      });
});

//var messages=["heatmaps of GPCR/G-protein","structure annotation of GPCR/G-protein"];
//var rank=0;
//
//// Code for Chrome, Safari and Opera
//document.getElementById("typewriter").addEventListener("webkitAnimationEnd", changeTxt);
//
//// Standard syntax
//document.getElementById("typewriter").addEventListener("animationend", changeTxt);
//
//function changeTxt(e){
//  _h1 = this.getElementsByTagName("h1")[0];
//  _h1.style.webkitAnimation = 'none'; // set element animation to none
//   setTimeout(function() { // you surely want a delay before the next message appears
//      _h1.innerHTML=messages[rank];
//      var speed = 3.5 * messages[rank].length/20; // adjust the speed (3.5 is the original speed, 20 is the original string length
//      _h1.style.webkitAnimation = 'typing '+speed+'s steps(40, end), blink-caret .75s step-end infinite'; //  switch to the original set of animation
//      (rank===messages.length-1)?rank=0:rank++; // if you have displayed the last message from the array, go back to the first one, else go to next message
//    }, 1000);
//}