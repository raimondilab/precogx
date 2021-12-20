$( document ).ready(function() {

/* Javascript to show and hide cookie banner using localstorage */
/* Shows the Cookie banner */
function showCookieBanner(){
 let cookieBanner = document.getElementById("cb-cookie-banner");
 cookieBanner.style.display = "block";
}

/* Hides the Cookie banner and saves the value to localstorage */
function hideCookieBanner(){
 localStorage.setItem("cb_isCookieAccepted", "yes");
 let cookieBanner = document.getElementById("cb-cookie-banner");
 cookieBanner.style.display = "none";
}

/* Checks the localstorage and shows Cookie banner based on it. */
function initializeCookieBanner(){
 let isCookieAccepted = localStorage.getItem("cb_isCookieAccepted");
 if(isCookieAccepted === null)
 {
  localStorage.setItem("cb_isCookieAccepted", "no");
  showCookieBanner();
 }
 if(isCookieAccepted === "no"){
  showCookieBanner();
 }
}

// Assigning values to window object
window.onload = initializeCookieBanner();
window.cb_hideCookieBanner = hideCookieBanner;


 $(document).on('click', '#ex1', function() {

      var seq1 = '# Dopamine Receptor 1\n';
      seq1 += ">DRD1\n";
      seq1 += "MRTLNTSAMDGTGLVVERDFSVRILTACFLSLLILSTLLGNTLVCAAVIRFRHLRSKVTN\n";
      seq1 += "FFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFGSFCNIWVAFDIMCSTASILNLCVISVD\n";
      seq1 += "RYWAISSPFRYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAKPTSPSDGNATSLA\n";
      seq1 += "ETIDNCDSSLSRTYAISSSVISFYIPVAIMIVTYTRIYRIAQKQIRRIAALERAAVHAKN\n";
      seq1 += "CQTTTGNGKPVECSQPESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCILPFCGS\n";
      seq1 += "GETQPFCIDSNTFDVFVWFGWANSSLNPIIYAFNADFRKAFSTLLGCYRLCPATNNAIET\n";
      seq1 += "VSINNNGAAMFSSHHEPRGSISKECNLVYLIPHAVGSSEDLKKEEAAGIARPLEKLSPAL\n";
      seq1 += "SVILDYDTDVSLEKIQPITQNGQHPT\n";
      seq1 += '# Free fatty acid receptor 2\n'
      seq1 += ">FFAR2\n";
      seq1 += "MLPDWKSSLILMAYIIIFLTGLPANLLALRAFVGRIRQPQPAPVHILLLSLTLADLLLLL\n";
      seq1 += "LLPFKIIEAASNFRWYLPKVVCALTSFGFYSSIYCSTWLLAGISIERYLGVAFPVQYKLS\n";
      seq1 += "RRPLYGVIAALVAWVMSFGHCTIVIIVQYLNTTEQVRSGNEITCYENFTDNQLDVVLPVR\n";
      seq1 += "LELCLVLFFIPMAVTIFCYWRFVWIMLSQPLVGAQRRRRAVGLAVVTLLNFLVCFGPYNV\n";
      seq1 += "SHLVGYHQRKSPWWRSIAVVFSSLNASLDPLLFYFSSSVVRRAFGRGLQVLRNQGSSLLG\n";
      seq1 += "RRGKDTAEGTNEDRGVGQGEGMPSSDFTTE";

      $('#validationTextarea').html(seq1);
    });

    $(document).on('click', '#ex2', function() {
      var seq2 = '# G-protein coupled receptor 183\n';
      seq2 += "P32249\n";
      seq2 += '# Thyrotropin receptor\n';
      seq2 += "TSHR_HUMAN";

      $('#validationTextarea').html(seq2);
    });

    $(document).on('click', '#ex3', function() {
      var seq3 = '# Free fatty acid receptor 2\n';
      seq3 += '#Complete loss of acetate-induced G protein-coupled receptor activity\n';
      seq3 += "FFAR2/Y90W/F61P\n";
      seq3 += "# G-protein coupled receptor 183\n";
      seq3 += "# 10-fold reduction in receptor activation (UniProt)\n";
      seq3 += "GPR183/Q287A";

      $('#validationTextarea').html(seq3);
    });

    $(document).on('click', '#ex4', function() {
      var seq4 = "# Parathyroid hormone/parathyroid hormone-related peptide receptor\n";
      seq4 += "PTH1R_HUMAN\n";
      seq4 += "# Metabotropic glutamate receptor 7 (Class C)\n";
      seq4 += "GRM7_HUMAN";

      $('#validationTextarea').html(seq4);
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
//}
