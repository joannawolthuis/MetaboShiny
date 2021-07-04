// https://colinfay.me/watch-r-shiny/
// http://usefulangle.com/post/45/animate-favicon-gif-javascript
// merged the two uwu

  var favicon_images = [
        'title1.png',
        'title2.png',
        'title3.png',
        'title4.png',
        'title5.png',
        'title6.png',
        'title7.png', 
        'title8.png',
        'title9.png'
        ];
        
        image_counter = 0; // To keep track of the current image

setInterval(function() {
    $("link[rel='icon']").remove();
    $("link[rel='shortcut icon']").remove();
    $("head").append('<link rel="icon" href="' + favicon_images[image_counter] + '" type="image/png">');
    var is_running = $('html').attr('class').includes('shiny-busy');
        if (is_running){
	// If last image then goto first image
	// Else go to next image    
	if(image_counter == favicon_images.length -1)
        image_counter = 0;
    else
        image_counter++;
          
        }else{
          $("link[rel='icon']").remove();
    $("link[rel='shortcut icon']").remove();
    $("head").append('<link rel="icon" href="' + favicon_images[0] + '" type="image/png">');
        }}, 700);
