$( document ).ready(function() {


// Assigning values to window object
window.onload = initializeCookieBanner();
window.hideCookieBanner = hideCookieBanner;


 $(document).on('click', '#ex1', function() {
      /*
      var seq1 = '# Dopamine Receptor 1\n';
      seq1 += ">DRD1\n";
      seq1 += "MRTLNTSAMDGTGLVVERDFSVRILTACFLSLLILSTLLGNTLVCAAVIRFRHLRSKVTN\n";
      seq1 += "FFVISLAVSDLLVAVLVMPWKAVAEIAGFWPFGSFCNIWVAFDIMCSTASILNLCVISVD\n";
      seq1 += "RYWAISSPFRYERKMTPKAAFILISVAWTLSVLISFIPVQLSWHKAKPTSPSDGNATSLA\n";
      seq1 += "ETIDNCDSSLSRTYAISSSVISFYIPVAIMIVTYTRIYRIAQKQIRRIAALERAAVHAKN\n";
      seq1 += "CQTTTGNGKPVECSQPESSFKMSFKRETKVLKTLSVIMGVFVCCWLPFFILNCILPFCGS\n";
      seq1 += "GETQPFCIDSNTFDVFVWFGWANSSLNPIIYAFNADFRKAFSTLLGCYRLCPATNNAIET\n";
      seq1 += "VSINNNGAAMFSSHHEPRGSISKECNLVYLIPHAVGSSEDLKKEEAAGIARPLEKLSPAL\n";
      seq1 += "SVILDYDTDVSLEKIQPITQNGQHPT\n\n";
      seq1 += '# Free fatty acid receptor 2\n'
      seq1 += ">FFAR2\n";
      seq1 += "MLPDWKSSLILMAYIIIFLTGLPANLLALRAFVGRIRQPQPAPVHILLLSLTLADLLLLL\n";
      seq1 += "LLPFKIIEAASNFRWYLPKVVCALTSFGFYSSIYCSTWLLAGISIERYLGVAFPVQYKLS\n";
      seq1 += "RRPLYGVIAALVAWVMSFGHCTIVIIVQYLNTTEQVRSGNEITCYENFTDNQLDVVLPVR\n";
      seq1 += "LELCLVLFFIPMAVTIFCYWRFVWIMLSQPLVGAQRRRRAVGLAVVTLLNFLVCFGPYNV\n";
      seq1 += "SHLVGYHQRKSPWWRSIAVVFSSLNASLDPLLFYFSSSVVRRAFGRGLQVLRNQGSSLLG\n";
      seq1 += "RRGKDTAEGTNEDRGVGQGEGMPSSDFTTE";

      $('#validationTextarea').html(seq1);
      */
      var seq1 = "# Example of alternative splicing in Thromboxane A2 receptor\n\n"
      seq1 += "# Canonical isoform of Thromboxane A2 receptor\n"
      seq1 += ">P21731|TA2R_HUMAN\n"
      seq1 += "MWPNGSSLGPCFRPTNITLEERRLIASPWFAASFCVVGLASNLLALSVLAGARQGGSHTR\n"
      seq1 += "SSFLTFLCGLVLTDFLGLLVTGTIVVSQHAALFEWHAVDPGCRLCRFMGVVMIFFGLSPL\n"
      seq1 += "LLGAAMASERYLGITRPFSRPAVASQRRAWATVGLVWAAALALGLLPLLGVGRYTVQYPG\n"
      seq1 += "SWCFLTLGAESGDVAFGLLFSMLGGLSVGLSFLLNTVSVATLCHVYHGQEAAQQRPRDSE\n"
      seq1 += "VEMMAQLLGIMVVASVCWLPLLVFIAQTVLRNPPAMSPAGQLSRTTEKELLIYLRVATWN\n"
      seq1 += "QILDPWVYILFRRAVLRRLQPRLSTRPRSLSLQPQLTQRSGLQ\n\n"
      seq1 += "# Isoform 2 of Thromboxane A2 receptor\n"
      seq1 += ">P21731-2|TA2R_HUMAN\n"
      seq1 += "MWPNGSSLGPCFRPTNITLEERRLIASPWFAASFCVVGLASNLLALSVLAGARQGGSHTR\n"
      seq1 += "SSFLTFLCGLVLTDFLGLLVTGTIVVSQHAALFEWHAVDPGCRLCRFMGVVMIFFGLSPL\n"
      seq1 += "LLGAAMASERYLGITRPFSRPAVASQRRAWATVGLVWAAALALGLLPLLGVGRYTVQYPG\n"
      seq1 += "SWCFLTLGAESGDVAFGLLFSMLGGLSVGLSFLLNTVSVATLCHVYHGQEAAQQRPRDSE\n"
      seq1 += "VEMMAQLLGIMVVASVCWLPLLVFIAQTVLRNPPAMSPAGQLSRTTEKELLIYLRVATWN\n"
      seq1 += "QILDPWVYILFRRAVLRRLQPRLSTRPRRSLTLWPSLEYSGTISAHCNLRLPGSSDSRAS\n"
      seq1 += "ASRAAGITGVSHCARPCMLFDPEFDLLAGVQLLPFEPPTGKALSRKD\n"

      $('#validationTextarea').html(seq1);
    });


    $(document).on('click', '#ex2', function() {
      var seq2 = '# G-protein coupled receptor 183\n';
      seq2 += "P32249\n\n";
      seq2 += '# Thyrotropin receptor\n';
      seq2 += "TSHR_HUMAN";

      $('#validationTextarea').html(seq2);
    });

    $(document).on('click', '#ex3', function() {
      var seq3 = '# Melanocyte-stimulating hormone receptor\n';
      seq3 += "# D294H mutation associated with a risk for developing melanoma;\n"
      seq3 += "# unable to stimulate cAMP production as strongly as the wild type receptor\n";
      seq3 += "# in response to alpha-melanocyte-stimulating hormone stimulation.\n";
      seq3 += "# PMID: 7581459, 11179997, 17616515, 17999355, 18366057, 19710684, 20585627\n";
      seq3 += "MC1R/D294H";

      $('#validationTextarea').html(seq3);
    });

    $(document).on('click', '#ex4', function() {
      var seq4 = "# Adhesion G protein-coupled receptor E5 (CD7)\n";
      seq4 += "# Known to couple GNA12 (Ward Y. et al, Cancer Res. 2011; PMID: 21978933)\n"
      seq4 += "ADGRE5";

      $('#validationTextarea').html(seq4);
    });

    $(document).on('click', '#ex5', function() {
      var seq5 = "# GtoPdb nomenclature\n";
      seq5 += '# Vasoactive intestinal peptide receptor 2\n';
      seq5 += "VPAC2 receptor";

      $('#validationTextarea').html(seq5);
    });

      //$('#validationTextarea').html("DRD1/F61P\nCXCR3/D278A");
      $("#running").hide();

      $("#submit").click(function(){

       var val = document.querySelector('#validationTextarea').value.trim();

       if (val != ""){
           $("#home").slideUp();

            var element = document.getElementById("footer");
            element.classList.add("footer");

            $("#running").show(1000);
        }
      });
});

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
