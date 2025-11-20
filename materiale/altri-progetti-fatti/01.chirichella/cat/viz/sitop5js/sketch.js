/**
 * Simulazione sistema + controllore in p5.js
 * Visibile su https://marziano.top/cat.html
 * con uno slider 0..360° per impostare il riferimento in deg.
 */

// -------------------------------------------------------------------
// 1) PARAMETRI
// -------------------------------------------------------------------
let k       = 500.0;          // [N/m^2]
let beta    = 0.5;            // [N*s/m^2]
let alphaDeg = 50;            // userà 50 gradi
let alpha = 0;                // verrà assegnato in setup() con radians(alphaDeg)
let J       = 400.0;          // [J*m^2]

// Equilibrio
let theta_eDeg = 0;           // in gradi
let theta_e = 0;              // in radianti
let Cm_e; 

// tau(theta)
function tauFunction(th) {
  // tau = cos(alpha) / [1 - (sin(alpha)*cos(theta))^2]
  let sa = sin(alpha);
  let ca = cos(alpha);
  return ca / (1.0 - (sa * cos(th)) * (sa * cos(th)));
}

//parametri dello slider
let sliderDeg = 36.0;  // Valore iniziale del riferimento in gradi
let sliderX = 10;
let sliderY = 130;
let sliderW = 200;  
let sliderH = 15;

// -------------------------------------------------------------------
// 2) Specifiche di controllo
// -------------------------------------------------------------------

// Vogliamo < e_star=0.01 for W=1.5, D_amp=1 => min_mu = 250
let R_s_val = 1.2290e+05;

console.log(R_s_val)

// R_d_antic(s) = (1 + tau_a s)/(1 + alpha_a tau_a s)
let tau_a   = 12.7498; 
let alpha_a = 2.2 * 10**-5;

// R_d_polo(s) = 1 / [ (1 + T1_f s) ]
let T1_f = 0.3 * 10^(-4);

// -------------------------------------------------------------------
// 3) Variabili di stato (Controller)
// -------------------------------------------------------------------

// A) R_d_antic: Y + alpha_a tau_a dY/dt = U + tau_a dU/dt nel tempo
let x_antic = 0.0;        
let input_antic_old = 0.0; // input precedente per calcolare derivata approx

// B) R_d_polo(s) = 1 / [(1 + T1_f s)]
let x_polo = 0.0;   // y(t)

// -------------------------------------------------------------------
// 4) Var di simulazione (Plant)
// -------------------------------------------------------------------
let theta = 0.0;
let omega = 0.0;

let dt    = 1e-7;  // Delta T
let t     = 0.0;
// riferimento in deg (slider), convertito in rad
let W;

// Disturbo
function disturbo(time) {
  // d(t) = 0.2*sin(0.1*1*t) + 0.2*sin(0.1*2*t) + 0.2*sin(0.1*3*t) + 0.2*sin(0.1*4*t);
  return 0.2 * sin(0.1 * 1.0 * time) 
       + 0.2 * sin(0.1 * 2.0 * time)
       + 0.2 * sin(0.1 * 3.0 * time)
       + 0.2 * sin(0.1 * 4.0 * time);
}

// -------------------------------------------------------------------
// 5) Parte grafica
// -------------------------------------------------------------------
let angleX = 0;  
let angleY = 0;          
let lastMouseX, lastMouseY;
let isDragging = false;  

let tableTex;

let motorAngle = 0.0;
let zoomLevel = 500;  // zoom

// Leggiamo lo slider dall'HTML
let angleSlider = null; 
let angleVal = null;  
let infoText = null;

function preload() {
  // Carica texture
  tableTex = loadImage("texture.png",
    () => { print("Texture caricata correttamente"); },
    () => { print("Errore nel caricare texture"); }
  );
}

function setup() {
  createCanvas(windowWidth, windowHeight , WEBGL);
  frameRate(60);

  // Trova l'elemento slider e il relativo label
  angleSlider = document.getElementById('angleSlider');
  angleVal = document.getElementById('angleVal');
  infoText = document.getElementById('infoText');

  // Imposta callback per cambiare sliderDeg
  angleSlider.addEventListener('input', () => {
    let val = angleSlider.value;
    sliderDeg = float(val);
    angleVal.innerText = val;
  });

  // Impostazioni iniziali
  alpha = radians(alphaDeg);
  theta_e = radians(theta_eDeg);

  Cm_e = (k * theta_e) / tauFunction(theta_e);


  theta = theta_e;
  omega = 0.0;

  x_antic = 0.0; 
  input_antic_old = 0.0; 
  x1_polo = 0.0; 
  x_polo = 0.0;

  angleX = radians(-20);
  angleY = radians(30);

  print("Tavolo rotante 3D - p5.js");
  print("Premi 'r' per resettare");
}

function draw() {
  background(220);

  noStroke();

  noLights()

  // Integra più volte per frame
  let stepsPerFrame = 100;
  for (let i = 0; i < stepsPerFrame; i++) {
    updateSystem();
  }

  push();
    
  translate(0, 0, zoomLevel);

  // Rotazioni della "camera"
  rotateX(angleX);
  rotateY(angleY);

  // Motore
  
  push();
  
  translate(-200, 100, 0);
  rotateZ(HALF_PI + radians(20));
  rotateY(motorAngle);

  drawCylinder(20, 30, 24, tableTex);
  
  stroke(0);
  strokeWeight(3);
  let topMotor = 15;
  let rodX = 0;
  let rodY = -213;
  line(0, topMotor, 0, rodX, rodY, 0);

  // Sfera giunto
  push();
  translate(rodX, rodY, 0);
  fill(255, 0, 0);
  noStroke();
  sphere(10);
  pop();
  pop();

  noStroke();
  push();
  rotateY(theta - theta_e);
  let tableRadius = 100;
  let tableHeight = 10;


  // Tavolo
  drawCylinder(tableRadius, tableHeight, 24, tableTex);
  pop();

  stroke(0);
  strokeWeight(3);
  line(rodX, rodY + 210, 0, rodX, rodY + 380, 0);
  pop();

  // Info testuali
  push();
  resetMatrix(); // torna in 2D
  fill(0);
  text("T = " + nf(t, 1, 5) + " s", 10, 20);
  text("theta (deg) = " + nf(degrees(theta % TWO_PI),1,2), 10, 40);
  text("Premi 'r' per resettare", 10, 60);
  pop();

  let infoString =
    "T = " + nf(t,1,5) + " s<br>" +
    "theta (deg) = " + nf(degrees(theta % TWO_PI),1,2) + "° <br>" +
    "Premi 'r' per resettare";
  infoText.innerHTML = infoString;

}

// Gestione zoom con rotellina mouse
function mouseWheel(event) {
  let e = event.delta;
  zoomLevel -= e * 0.2;  
}

// Aggiorno il sistema con il metodo di Eulero
function updateSystem() {
  // Convert sliderDeg => rad
  W = radians(sliderDeg);

  //segnale di errore
  let e = W - theta;
  //regolazione statica
  let e_dist = R_s_val * e;

  //input del reg dinamico
  let u_antic = e_dist;
  //derivata dell'input
  let du_antic = (u_antic - input_antic_old) / dt;

  //aggiorno storico
  input_antic_old = u_antic;

  //derivata discreta nel tempo dell'anticipatore
  let x_antic_dot = (u_antic + tau_a * du_antic - x_antic) / (alpha_a * tau_a);
  x_antic += x_antic_dot * dt;
  //output anticipato, input del polo
  let out_antic = x_antic;
  let in_polo = out_antic;

  //derivata nel tempo del polo
  let polo_dot = (in_polo - x_polo)  / (T1_f);
  //aggiorno stato post polo
  x_polo += polo_dot * dt;
  let out_polo = x_polo;
  let coppia = Cm_e + out_polo;


  //let dist = disturbo(t);
  let dist = 0;
  //sistema non lineare
  let dtheta = omega;
  let domega = (tauFunction(theta)*coppia - beta*omega - k*theta + dist) / J;
  theta += dtheta * dt;
  omega += domega * dt;

  //simulo la curva rpm coppia
  motorAngle += 1e-5 * coppia * dt;
  t += dt;
}

// Disegna un cilindro
function drawCylinder(radius, h, segments, topTexture) {
  let halfH = h / 2;

  // top face
  if (topTexture) {
    textureMode(NORMAL);
    beginShape(TRIANGLE_FAN);
    texture(topTexture);
    vertex(0, halfH, 0, 0.5, 0.5);
    for (let i = 0; i <= segments; i++) {
      let angle = TWO_PI * i / segments;
      let x = radius * cos(angle);
      let z = radius * sin(angle);
      let u = 0.5 + 0.5 * cos(angle);
      let v = 0.5 + 0.5 * sin(angle);
      vertex(x, halfH, z, u, v);
    }
    endShape();
  } else {
    fill(180, 180, 220);
    beginShape(TRIANGLE_FAN);
    vertex(0, halfH, 0);
    for (let i = 0; i <= segments; i++) {
      let angle = TWO_PI * i / segments;
      let x = radius * cos(angle);
      let z = radius * sin(angle);
      vertex(x, halfH, z);
    }
    endShape();
  }

  // bottom face
  fill(180,180,220);
  beginShape(TRIANGLE_FAN);
  vertex(0, -halfH, 0);
  for (let i = segments; i >= 0; i--) {
    let angle = TWO_PI * i / segments;
    let x = radius * cos(angle);
    let z = radius * sin(angle);
    vertex(x, -halfH, z);
  }
  endShape();

  if (topTexture) {
    beginShape(TRIANGLE_FAN);
    texture(topTexture);
    textureMode(NORMAL);
    vertex(0, -halfH, 0, 0.5, 0.5);
    for (let i = segments; i >= 0; i--) {
      let angle = TWO_PI * i / segments;
      let x = radius * cos(angle);
      let z = radius * sin(angle);
      let u = 0.5 + 0.5 * cos(angle);
      let v = 0.5 + 0.5 * sin(angle);
      vertex(x, -halfH, z, u, v);
    }
    endShape();
  }

  // lato
  fill(180,180,220);
  beginShape(QUAD_STRIP);
  for (let i = 0; i <= segments; i++) {
    let angle = TWO_PI * i / segments;
    let x = radius * cos(angle);
    let z = radius * sin(angle);
    // top vertex
    vertex(x, halfH, z);
    // bottom vertex
    vertex(x, -halfH, z);
  }
  endShape();
}

// Gestione slider con mouse
function mousePressed() {
  let sx = sliderX, sy = sliderY, sw = sliderW, sh = sliderH;

  if (
    mouseX >= sx && mouseX <= sx + sw &&
    mouseY >= sy && mouseY <= sy + sh
  ) {
  } else {
    isDragging = true;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
  }
}

function mouseReleased() {
  isDragging = false;
}

function mouseDragged() {
  if (!isDragging) return;
  let dx = mouseX - lastMouseX;
  let dy = mouseY - lastMouseY;
  angleY += dx * 0.01;
  angleX -= dy * 0.01;
  lastMouseX = mouseX;
  lastMouseY = mouseY;
}

// Handler per 'r'
function keyPressed() {
  if (key === 'r' || key === 'R') {
    t = 0.0;
    theta = radians(theta_eDeg);
    omega = 0.0;
    x_antic = 0.0; 
    input_antic_old = 0.0; 
    x1_polo = 0.0; 
    x_polo = 0.0;
    angleX = radians(-20);
    angleY = radians(30);
    sliderDeg = 36.0;
    if (angleSlider) {
      angleSlider.value = 36;
      document.getElementById('angleVal').innerText = '36';
    }
    motorAngle = 0.0;
  }
}

