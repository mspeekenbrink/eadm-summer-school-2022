remark.macros.scale_img = function (percentage) {
  var url = this;
  return '<img src="' + url + '" style="width: ' + percentage + '" />';
};

remark.macros.scale_text = function (percentage) {
  var text = this;
  return '<div style="font-size: ' + percentage + '; " />' + text + '</div>';
};

remark.macros.image_source = function () {
  var url = this;
  return '<div style="font-size: .6rem; opacity: 0.8">Image source: ' + url + '</div>';
};

remark.macros.image_source_footer = function () {
  var url = this;
  return '<div style="font-size: .6rem; opacity: 0.8; position: abolute; bottom:0">Image source: ' + url + '</div>';
};

remark.macros.position_image = function(width="30%", left="85%", top="15%", opacity="1", cl="") {
var url = this;
return '<img src="' + url + '" style="position:absolute;left:' + left + ';top:' + top + ';width:' + width + ';opacity:' + opacity + '" class="' + cl + '" />';
};

remark.macros.position_text = function(width="30%", left="85%", top="15%", cl="") {
var text = this;
return '<div style="position:absolute;left:' + left + ';top:' + top + ';width:' + width + '" class="' + cl + '" />' + text + '</div>';
};