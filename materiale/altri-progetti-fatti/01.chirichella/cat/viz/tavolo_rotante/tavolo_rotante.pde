/**
 * Simulazione sistema + controllore in Processing.js
 * con uno slider 0..360° per impostare il riferimento in deg.
 */

// -------------------------------------------------------------------
// 1) PARAMETRI
// -------------------------------------------------------------------
float k       = 500.0;          // [N/m^2]
float beta    = 0.5;            // [N*s/m^2]
float alpha   = radians(50.0);  // [rad]
float J       = 400.0;          // [J*m^2]

// Equilibrio
float theta_e = radians(0);
float Cm_e; 

// tau(theta)
float tauFunction(float th) {
  // tau = cos(alpha) / [1 - (sin(alpha)*cos(theta))^2]
  float sa = sin(alpha);
  float ca = cos(alpha);
  return ca / (1.0 - (sa * cos(th)) * (sa * cos(th)));
}

// -------------------------------------------------------------------
// 2) Specifiche di controllo
// -------------------------------------------------------------------

// Vogliamo < e_star=0.01 for W=1.5, D_amp=1 => min_mu = 250
float min_mu  = 250.0;
// R_s_val = 21.5 * min_mu => 21.5 * 250 = 5375
float R_s_val = 21.5 * min_mu;

// R_d_dist = 50000
float R_d_dist = 50000.0;

// R_d_antic(s) = (1 + tau_a s)/(1 + alpha_a tau_a s)
float tau_a   = 0.5e-2;   // 0.005
float alpha_a = 0.01;

// R_d_polo(s) = 1 / [ (1 + T1_f s)(1 + T2_f s ) ]
float T1_f = 1.0e-4;
float T2_f = 0.8e-4;

// -------------------------------------------------------------------
// 3) Variabili di stato (Controller)
// -------------------------------------------------------------------

// A) R_d_antic: Y + alpha_a tau_a dY/dt = U + tau_a dU/dt
float x_antic = 0.0;        
float input_antic_old = 0.0; // input precedente per calcolare derivata approx

// B) R_d_polo(s) = 1 / [(1 + T1_f s)(1 + T2_f s)]
float x1_polo = 0.0;   // y(t)
float x2_polo = 0.0;   // y'(t)

// -------------------------------------------------------------------
// 4) Var di simulazione (Plant)
// -------------------------------------------------------------------
float theta = 0.0;
float omega = 0.0;

float dt    = 1e-7;  // Delta T
float t     = 0.0;
//ref
float W;

// Disturbo
float disturbance(float time) {
  // d(t) = 0.2*sin(0.1*1*t) + 0.2*sin(0.1*2*t) + ...
  return 0.2 * sin(0.1f * 1.0f * time) 
       + 0.2 * sin(0.1f * 2.0f * time)
       + 0.2 * sin(0.1f * 3.0f * time)
       + 0.2 * sin(0.1f * 4.0f * time);
}

// -------------------------------------------------------------------
// 5) Parte grafica
// -------------------------------------------------------------------
float angleX = radians(-20);  
float angleY = radians(30);          
float lastMouseX, lastMouseY;
boolean isDragging = false;  

PImage tableTex;

float sliderDeg = 360.0; 

int sliderX = 10;
int sliderY = 130;
int sliderW = 200;  
int sliderH = 15;

boolean sliderDragging = false;

float motorAngle = 0.0;

// -------------------------------------------------------------------
// 6) Setup
// -------------------------------------------------------------------
void setup() {
  size(1000, 600, P3D);
  frameRate(60);

  tableTex = loadImage("texture.png");
  if (tableTex == null) {
    println("Errore nel caricare texture");
  }

  Cm_e = (k * theta_e) / tauFunction(theta_e);

  theta = theta_e;
  omega = 0.0;
  
  x_antic = 0.0; 
  input_antic_old = 0.0; 
  x1_polo = 0.0; 
  x2_polo = 0.0;

  println("Tavolo rotante 3D");
  println("Premo 'r' per resettare e cambia il ref dallo slider");
}


float zoomLevel = 100;  

// -------------------------------------------------------------------
// 8) loop grafico
// -------------------------------------------------------------------
void draw() {
  background(255);
  lights();
  noStroke();

  int stepsPerFrame = 100;  
  for (int i = 0; i < stepsPerFrame; i++) {
    updateSystem();
  }

  pushMatrix();
  translate(width / 2, height / 2, zoomLevel);  
  rotateX(angleX);
  rotateY(angleY);

  // Motore
  pushMatrix();
  translate(-200, 100, 0);
  rotateZ(HALF_PI + radians(20));
  rotateY(motorAngle);
  emissive(50, 50, 50);
  drawCylinder(20, 30, 24, tableTex);

  stroke(0);
  strokeWeight(3);
  float topMotor = 15;
  float rodX = 0;
  float rodY = -213;
  line(0, topMotor, 0, rodX, rodY, 0);

  pushMatrix();
  translate(rodX, rodY, 0);
  fill(255, 0, 0);
  noStroke();
  sphere(10);
  popMatrix();
  popMatrix();

  noStroke();
  pushMatrix();
  rotateY(theta - theta_e);
  float tableRadius = 100;
  float tableHeight = 10;
  emissive(50, 50, 50);
  drawCylinder(tableRadius, tableHeight, 24, tableTex);
  popMatrix();

  stroke(0);
  strokeWeight(3);
  line(rodX, rodY + 210, 0, rodX, rodY + 380, 0);
  popMatrix();

  fill(0);
  text("T = " + nf(t,1,5) + " s", 10, 20);
  text("theta (deg) = " + nf(degrees(theta % TWO_PI),1,2), 10, 40);
  text("Premi 'r' per resettare", 10, 60);
  text("Usa lo slider per impostare ref (deg)", 10, 80);

  drawSlider();
}

// -------------------------------------------------------------------
// Handler per mouse wheel (zoom)
// -------------------------------------------------------------------
void mouseWheel(MouseEvent event) {
  float e = event.getCount();
  zoomLevel -= e * 20;  
}


// -------------------------------------------------------------------
// Aggiorno il sistema con il metodo di Eulero
// -------------------------------------------------------------------
void updateSystem() {
  W = radians(sliderDeg);
  float e = W - theta;
  float e_s = R_s_val * e;
  float e_dist = R_d_dist * e_s;
  float u_antic = e_dist;
  float du_antic = (u_antic - input_antic_old) / dt;
  input_antic_old = u_antic;
  float x_antic_dot = (u_antic + tau_a*du_antic - x_antic) / (alpha_a * tau_a);
  x_antic += x_antic_dot * dt;
  float out_antic = x_antic;
  float in_polo = out_antic;
  float y_pp = (in_polo - x1_polo - (T1_f + T2_f)*x2_polo) / (T1_f * T2_f);
  x2_polo += y_pp * dt;
  x1_polo += x2_polo * dt;
  float out_polo = x1_polo;
  float coppia = Cm_e + out_polo;
  float dist = disturbance(t);
  float dtheta = omega;
  float domega = (tauFunction(theta)*coppia - beta*omega - k*theta + dist) / J;
  theta += dtheta * dt;
  omega += domega * dt;
  
  //rallento il motore per una visualizzazione migliore : posso immaginare che ci siano degli ingranaggi tra il motore e il giunto
  motorAngle += 1e-5 * coppia * dt;
  t += dt;
}

// -------------------------------------------------------------------
// grafica per disegnare tavolo
// -------------------------------------------------------------------
void drawCylinder(float radius, float h, int segments, PImage topTexture) {
  float halfH = h/2;
  if (topTexture != null) {
    textureMode(NORMAL);
    beginShape(TRIANGLE_FAN);
    texture(topTexture);
    vertex(0, halfH, 0, 0.5, 0.5);
    for(int i=0; i<=segments; i++){
      float angle = TWO_PI * i / segments;
      float x = radius * cos(angle);
      float z = radius * sin(angle);
      float u = 0.5 + 0.5 * cos(angle);
      float v = 0.5 + 0.5 * sin(angle);
      vertex(x, halfH, z, u, v);
    }
    endShape();
  }
  fill(180,180,220);
  beginShape(TRIANGLE_FAN);
  vertex(0, -halfH, 0); 
  for(int i=segments; i>=0; i--){
    float angle = TWO_PI * i / segments;
    float x = radius * cos(angle);
    float z = radius * sin(angle);
    vertex(x, -halfH, z);
  }
  endShape();
  beginShape(TRIANGLE_FAN);
  if (topTexture != null) {
    texture(topTexture);
    textureMode(NORMAL);
    vertex(0, -halfH, 0, 0.5, 0.5);
    for(int i=segments; i>=0; i--){
      float angle = TWO_PI * i / segments;
      float x = radius * cos(angle);
      float z = radius * sin(angle);
      float u = 0.5 + 0.5 * cos(angle);
      float v = 0.5 + 0.5 * sin(angle);
      vertex(x, -halfH, z, u, v);
    }
  }
  endShape();
  beginShape(QUAD_STRIP);
  for(int i=0; i<=segments; i++){
    float angle = TWO_PI * i / segments;
    float x = radius * cos(angle);
    float z = radius * sin(angle);
    vertex(x,  halfH, z);
    vertex(x, -halfH, z);
  }
  endShape();
}

// -------------------------------------------------------------------
// Codice per lo slider
// -------------------------------------------------------------------
void drawSlider() {
  stroke(0);
  fill(220);
  rect(sliderX, sliderY, sliderW, sliderH);
  float handleX = map(sliderDeg, 0, 360, 0, sliderW);
  fill(100, 100, 250);
  rect(sliderX + handleX - 5, sliderY - 5, 10, sliderH + 10);
  fill(0);
  text("Ref (deg) = " + nf(sliderDeg,1,1), sliderX + sliderW + 10, sliderY + 12);
}

// -------------------------------------------------------------------
// Handler per 'r'
// -------------------------------------------------------------------
void keyPressed() {
  if (key == 'r') {
    t = 0.0;
    theta = theta_e;
    omega = 0.0;
    x_antic = 0.0; 
    input_antic_old = 0.0; 
    x1_polo = 0.0; 
    x2_polo = 0.0;
    angleX = radians(-45);
    angleY = 0.0;
    sliderDeg = 360.0;
    motorAngle = 0.0;
  }
}

// -------------------------------------------------------------------
// Input mouse
// -------------------------------------------------------------------
void mousePressed() {
  if (mouseX >= sliderX && mouseX <= sliderX + sliderW &&
      mouseY >= sliderY && mouseY <= sliderY + sliderH) {
    sliderDragging = true;
    updateSliderValue(mouseX);
  } else {
    isDragging = true;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
  }
}

void mouseReleased() {
  isDragging = false;
  sliderDragging = false;
}

void mouseDragged() {
  if (sliderDragging) {
    updateSliderValue(mouseX);
  }
  else if (isDragging) {
    float dx = mouseX - lastMouseX;
    float dy = mouseY - lastMouseY;
    angleY += dx * 0.01;
    angleX -= dy * 0.01;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
  }
}

void updateSliderValue(float mx) {
  float clamped = constrain(mx, sliderX, sliderX + sliderW);
  sliderDeg = map(clamped, sliderX, sliderX + sliderW, 0, 360);
}
