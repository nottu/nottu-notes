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
                page.viewportSize = { width: 1920, height: 1080 };

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
                  var name = "capture_" + page.url + i + ".png";
                  console.log(name);
                  page.render(name);
                }
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
var start_url = 'https://getbootstrap.com/';
new phantom.Crawler().crawl(start_url, 2, 
    function onSuccess(page) {
        console.log("Loaded page. URL = " + page.url + " content length = " + page.content.length + " status = " + page.status);
    }, 
    function onFailure(page) {
        console.log("Could not load page. URL = " +  page.url + " status = " + page.status);
    }
);