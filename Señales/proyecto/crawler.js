(function(host) {

    function Crawler() {
        this.visitedURLs = {};
    };
    
    Crawler.webpage = require('webpage');

    Crawler.prototype.crawl = function (url, depth, onSuccess, onFailure) {
        if (0 == depth || this.visitedURLs[url]) {
            return;
        };
        var self = this;
        var page = Crawler.webpage.create();
        // console.log(typeof url)
        if( url.indexOf('github') != -1 ){
            console('nope')
            return;
        }
        page.open(url, function (status) {
            if ('fail' === status) { 
                onFailure({
                    url: url, 
                    status: status
                });
            } else {
                var documentHTML = page.evaluate(function () {
                    return document.body && document.body.innerHTML ? document.body.innerHTML : "";
                });
                self.crawlURLs(self.getAllURLs(page), depth - 1, onSuccess, onFailure);
                self.visitedURLs[url] = true;
                onSuccess({
                    url: url,
                    status: status,
                    content: documentHTML
                });
            };
        });
    };

    Crawler.prototype.getAllURLs = function(page) {
        return page.evaluate(function () {
            return Array.prototype.slice.call(document.querySelectorAll("a"), 0)
                .map(function (link) {
                    return link.getAttribute("href");
                });
        });
    };

    Crawler.prototype.crawlURLs = function(urls, depth, onSuccess, onFailure) {
        var self = this;
        urls.filter(function (url) {
            return Crawler.isTopLevelURL(url);
        }).forEach(function (url) {
            self.crawl(url, depth, onSuccess, onFailure);
        });
    };

    Crawler.isTopLevelURL = function(url) {
        return 0 == url.indexOf("http");
    };

    host.Crawler = Crawler;
})(phantom);


// var start_url = 'https://materializecss.com/';
// start_url = 'https://material-ui.com/'
// start_url = 'https://getbootstrap.com/'
// start_url = 'https://material-ui.com/demos/buttons/'
// start_url = 'https://semantic-ui.com/'
start_url = 'https://codepen.io/'

num = 50;
new phantom.Crawler().crawl(start_url, 2, 
    function onSuccess(pg) {
        num += 1
        console.log("Loaded page. URL = " + pg.url + " content length = " + pg.content.length + " status = " + pg.status);

        var page = require('webpage').create();
        var margin = 10;
        var selector = 'button, .btn';

        prename = 'button_' + num;

        page.open(pg.url, function(status) {
          page.viewportSize = { width: 1920, height: 1080 };
          console.log("Status: " + status);
          if(status === "success") {
            var rects = page.evaluate(function(){
              var items = document.querySelectorAll('.button, .btn, .st-btn, .c-btn');
              var rects = []
              for (var i = items.length - 1; i >= 0; i--) {
                var clipRect = items[i].getBoundingClientRect();
                if(clipRect.height > 0 && clipRect.width > 0) rects.push(clipRect);
              }
              return rects;
            });
            // var rects = getClipingRects(this, selectors)
            console.log(rects.length, 'items');
            for (var i = rects.length - 1; i >= 0; i--) {
              var clipRect = rects[i];
              page.clipRect = {
                  top:    clipRect.top - margin,
                  left:   clipRect.left - margin,
                  width:  clipRect.width + 2 * margin,
                  height: clipRect.height + 2 * margin
              };
              var name =  prename + "_" + i + ".png";
              console.log(name);
              page.render(name);
            }
          }
        });


    }, 
    function onFailure(page) {
        console.log("Could not load page. URL = " +  page.url + " status = " + page.status);
    }
);
