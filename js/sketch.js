// Main Sketch
let screensize = (4 * $("#jumbo-canvas").width()) / 5;
let scale;
let shdr;
let zscl = 100000.0;
//------------------------------------------------------------------------------
preload = function()
{
  shdr = loadShader('/js/shaders/vert.shader', '/js/shaders/frag.shader');
};
//------------------------------------------------------------------------------
function setup()
{
  var canvas = createCanvas(screensize, screensize, WEBGL);
  frameRate(30);
  // gl = canvas.getContext('webgl');
  canvas.parent('sketch-holder');
  pg = createGraphics(200, 200);
  pg.textSize(75);
  pg.background(0, 100);
  fill(255);
  texture(pg);
  shader(shdr);
  drawMesh();
  // console.log(portouts[40][40][1]);
}
//------------------------------------------------------------------------------
function draw()
{
  drawMesh();

  if (frameCount % 30 === 0)
  {
      console.log(portouts[40][40][1]);
  }
}
//------------------------------------------------------------------------------
function windowResized()
{
  screensize = (4 * $("#jumbo-canvas").width()) / 5;
  resizeCanvas(screensize, screensize);
}
//------------------------------------------------------------------------------
function drawMesh()
{
  updateState();

  background(0);
  scale = width / (Njx - 1);
  push();
  rotateX(PI / 3);
  rotateZ(millis() * 0.0002);
  rotateY(millis() * 0.0001);
  translate(-width / 2, -height / 2, 0);
  for (var yi = 0; yi < Njy - 1; ++yi)
  {
    beginShape(TRIANGLE_STRIP);
    for (var xi = 0; xi < Njx; ++xi)
    {
      vertex(xi * scale, yi * scale, Math.floor(portouts[xi][yi][1] * zscl));
      vertex(xi * scale, (yi + 1) * scale, Math.floor(portouts[xi][yi + 1][1] * zscl));
    }
    endShape()
  }
  pop();
}
//------------------------------------------------------------------------------
function mousePressed()
{}
