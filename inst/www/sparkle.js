(function($, window, document) {

  const defaults = {
    fill: "#fff",
    stroke: "#000",
    size: 20,
    delay: 0,
    duration: 1500,
    pause: 1000
  };
  const optionsKeys = ["fill", "stroke", "size", "delay", "duration", "pause"];
  
  const optionsStrToObj = function(optionsStr) {
    const optionsArr = !!optionsStr ? optionsStr.split(" ") : [];
    var optionsObj = {};
    
    for (var i=0; i<optionsArr.length; ++i) {
      optionsObj[optionsKeys[i]] = optionsArr[i];
    }
    
    return optionsObj;
  };
  
  $.fn.sparkle = function(options) {
    $.destroySparkle = $.destroySparkle || {};
    
    const $this = this;
    const id = this.data("sparkle-id") || (new Date()).getTime() + Math.random();
    
    if (options === "destroy" && this.find("svg").length > 0) {
      $.destroySparkle[id] = true;
      this.data("sparkle-id", null);
    }
    
    var settings;

    if (options instanceof Array) {
      for (var i=0; i<options.length; ++i) {
        $this.sparkle(optionsStrToObj(options[i]));
      }
      
      return;
    } else if (options instanceof Object) {
      settings = $.extend({}, defaults, options);
    } else {
      settings = defaults;
    }
    
    const cssAnimationAttr = "my-sparkle " + settings.duration + "ms infinite linear";

    const $star = $('<svg class="my-sparkle" version="1.1" viewBox="0.0 0.0 50.0 50.0" fill="none" stroke="none" stroke-linecap="square" stroke-miterlimit="10" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"><clipPath id="p.0"><path d="m0 0l50.0 0l0 50.0l-50.0 0l0 -50.0z" clip-rule="nonzero"></path></clipPath><g clip-path="url(#p.0)"><path fill="' + settings.stroke + '" fill-opacity="0.0" d="m0 0l50.0 0l0 50.0l-50.0 0z" fill-rule="nonzero"></path><path fill="' + settings.fill + '" d="m0.62204725 25.0l20.068499 -4.323374l4.309454 -20.13332l4.309454 20.13332l20.068499 4.323374l-20.068499 4.323374l-4.309454 20.133318l-4.309454 -20.133318z" fill-rule="nonzero"></path><path stroke="' + settings.stroke + '" stroke-width="1.0" stroke-linejoin="round" stroke-linecap="butt" d="m0.62204725 25.0l20.068499 -4.323374l4.309454 -20.13332l4.309454 20.13332l20.068499 4.323374l-20.068499 4.323374l-4.309454 20.133318l-4.309454 -20.133318z" fill-rule="nonzero"></path></g></svg>').css({
        position: "absolute",
        width: settings.size,
        height: settings.size,
        zIndex: 9999
    });

    /*const w = $('.sparkley').width();
    const h = $('.sparkley').height();*/
    w = 300;
    h = 25;

    const getCoordinates = function() {
      return {
        left: Math.random() * w,
        top: Math.random() * h
      };
    };

    const placeStar = function(init) {
      const coords = getCoordinates();

      if (init) {
        $this.append($star);
      }

      $star.css({
        "-moz-animation": cssAnimationAttr,
        "-webkit-animation": cssAnimationAttr,
        animation: cssAnimationAttr,
        display: "block",
        left: coords.left,
        top: coords.top
      });

      window.setTimeout(function() {
        $star.css({
          "-moz-animation": null,
          "-webkit-animation": null,
          animation: null,
          display: "none"
        });
        
        if (!$.destroySparkle[id]) {
          window.setTimeout(function() {
            placeStar(false);
          }, settings.pause);
        } else {
          $star.remove();
        }
      }, settings.duration);
    };

    if (this.css("position") === "static") {
      this.css("position", "relative");
    }

    if (!$.destroySparkle[id] && options !== "destroy") {
      window.setTimeout(function() {
        placeStar(true);
      }, settings.delay);

      this.data("sparkle-id", id);
    }

    return this;
  };
  
  $("#sparkley").sparkle({
    size: 25,
  }).sparkle({
    delay: 1000,
    pause: 750,
    size: 10
  });

  /*window.setTimeout(function() {
    $("#sparkley").sparkle("destroy");
  }, 21000);*/

})(jQuery, window, document);