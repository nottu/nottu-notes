var page = require('webpage').create();
var url = 'https://getbootstrap.com/docs/4.1/components/buttons/';
var margin = 10;
var selectors = ['button']


page.open(url, function(status) {
  page.viewportSize = { width: 1920, height: 1080 };
  console.log("Status: " + status);
  if(status === "success") {
    console.log('h1');
    var rects = page.evaluate(function(){
      var items = document.querySelectorAll('button, .btn');
      console.log(items.length, 'items')
      var rects = []
      for (var i = items.length - 1; i >= 0; i--) {
        var clipRect = items[i].getBoundingClientRect();
        if(clipRect.height > 0 && clipRect.width > 0) rects.push(clipRect);
      }
      return rects;
    });
    // var rects = getClipingRects(this, selectors)
    console.log(rects.length);
    for (var i = rects.length - 1; i >= 0; i--) {
      var clipRect = rects[i];
      page.clipRect = {
          top:    clipRect.top - margin,
          left:   clipRect.left - margin,
          width:  clipRect.width + 2 * margin,
          height: clipRect.height + 2 * margin
      };
      console.log(page.clipRect.height)
      var name = "capture_" + i + ".png";
      console.log(name);
      page.render(name);
      console.log('h5');
    }
  }
  phantom.exit();
});