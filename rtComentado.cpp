// para cambiar de muestreo es con las funciones shade
//  ************************** DOMINGO *****************************************************
#include <math.h> // se usa la función sqrt(), cos(), sin(), acos(), asin(), exp(), log(), pow(), fabs(), fmin(), fmax(), y M_PI
#include <stdlib.h> // se usa rand(), RAND_MAX, srand(). Principalmente se usa para los números aleatorios
#include <stdio.h>  
#include <omp.h> // permite la parelelización
#include <time.h> // se usa time(). Semilla para los números aleatorios


#define N_ 32 // es el número de muestras. Se suelen usar en potencias de 2 
#define N_ESFERAS 9 // el número de esferas en la escena
#define N_ARREGLO 8 // el tamaño del arreglo

class Vector // se declara la clase vector
{
public: // todos los parámetros son públicos
	double x, y, z; // coordenadas x,y,z 

    /*
    La suma se usa para componer direcciones, sumar contribuciones de luz, mover puntos
    La resta para calcular vectores de dirección
    El * para escalar vecctores, para aplicar las BRDF
    El producto cruz se usa calcular normales, para hacer vectores ortogonales, como en CoordinateSystem
    El producto punto es el escalar entre dos vectores. Se usa para cálcular el coseno de un ángulo, la intensidad de lu< en una superficie, y cálculos de intersección
    Mult se usa en montecarlo y montecarlofuentepuntual. Es para modular colore sde luz 
    normalice normaliza vectores jsjs y se usa en muchas funiones, como localAGlobal(), coodinateSsytem(), intersecc()
    */
  
	// Constructor del vector, parametros por default en cero
	Vector(double x_= 0, double y_= 0, double z_= 0){ x=x_; y=y_; z=z_; }
  
	// operador para suma y resta de vectores
	Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
	Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
	// operator multiplicacion vector y escalar 
	Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }
  
	// operator % para producto cruz
	Vector operator%(Vector&b){return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
	
	// producto punto con vector b
	double dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; }

	// producto elemento a elemento (Hadamard product)
	Vector mult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
	
	// normalizar vector 
	Vector& normalize(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
};
typedef Vector Point; // es un alias ya que se puede hacer lo mismo con un vector, pero viendolo desde un punto distinto. 
typedef Vector Color; // lo mismo para color, es un alias ya que puede hacer lo mismo que vector pero desde otro punto de vista

enum MaterialType { // definimos una enumeración con los tipos de materiales que tenemos
    DIFFUSO, // 0
    CONDUCTOR // 1
};

class Ray // definimos la clase rayo
{ 
public: // todos los atributos van a ser púplicos
	Point o; // el punto de origen
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};
 
class Sphere  // definimos la clase esfera
{
public: // los tributos van a ser públicos
	double r;	// radio de la esfera
	Point p;	// posicion central de la esfera
	Color c;	// color  
	Color e;  // la emisión de la esfera
    MaterialType material; //el tipo de material, difuso o conductor
    double aspereza; // para microfacet 
    Color eta, kappa; // para conductor. ETA índice de refracción y dice como el material ralentiza la luz, para cálcular como la luz se refleja.   
                      //   KAPPA coeficiente de extinción, describe como el material absorve la luz, cálcula la reflectancia.
	Sphere(double r_, Point p_, Color c_, Color e_ = Color(), MaterialType mat = DIFFUSO, double aspe = 0.1, Color eta_ = Color(1,1,1), Color kappa_ = Color(0,0,0)) 
        : r(r_), p(p_), c(c_), e(e_), material(mat), aspereza(aspe), eta(eta_), kappa(kappa_) {}
   
	double intersect(const Ray &ray) const { // se define la función intersect en 80 y recibe un rayo constante
		Vector op = p - ray.o; // del origen del rayo al centro de la esfera
		Vector d = ray.d; // dirección del rayo
		double eps = 1e-4; // número chikito para evitar intersecciones fantasma
		double b = op.dot(d); // es es producto punto del rayo entre el origen y el centro de la esfera por la dirección del rayo
		double c = op.dot(op) - r*r; // la otra parte de la ecuación
		double delta = b*b - c; // para ver si hay intersección
		if (delta < 0) { return 0.0; } // si es menor, no hay intersección
		delta = sqrt(delta); 

		double t = b - delta;
		if (t > eps) { return t; } // se evalúa primero la resta por que sería la más cercana
		t = b + delta; // si fue negariva, se evalúa la suma, que sería t2
		return (t >  eps) ? t : 0.0; // si sale negativa significa que la intersección sucedió atrás del rayo, por lo que no es válido
	}
};
// --------------------------- HORA 1 ------------------------------------------------------------------------
class FuenteLuminosa : public Sphere { // se define la clase de fuente luminosa, que hereda de esfera en 97
public: // susatributos son públicos
	FuenteLuminosa(double r_, Point p_, Color c_, Color e_ = Color(10, 10,10)) //la diferencia es que le ponemos que la emisión sea fina de 10, auqne se podría cambiar al crear un objeto
		: Sphere(r_, p_, c_, e_) {}	// se definen sus valores
};

class FuentePuntual : public Sphere { // se define otra clase heredada de esfera en 103. En ella lo que se cambia es la emisión y el radio
public: // los atributos son públicos
    FuentePuntual(Point p_, Color e_ = Color(4000, 4000, 4000), double r_ = 0)  // la emisión debe de ser muy grande ya es que sólo un punto el que da luz a toda la escena
        : Sphere(r_, p_, Color(), e_) {} 
};

struct MuestreoResult { // se crea una estructura para manejar el resultado de los muestreos, ya que de todos se ocupa casi lo mismo
	Vector wi; // un vector entrante
	double prob; // la probabilidad de la intersección
	Point puntoEnFuente; // el punto donde toca una fuente de luz. sólo se termina usando en muestreo de área, fuente puntual y ángulo sólido

	MuestreoResult(Vector w, double p, Point pt = Point()) : wi(w), prob(p), puntoEnFuente(pt) {} // el punto en fuente por defecto es 0
};

Sphere spheres[] = {
    //Escena: radio, posicion, color, emisión
    // para las esferas normales, por defecto la emisión es 0, y la aspereza 1,eta 1 y kappa 0. estos son los valores para materiales difusos por defecto
     	Sphere(1e5,  Point(-1e5 - 49, 0, 0),     Color(.75, .25, .25)), // pared izq
        Sphere(1e5,  Point(1e5 + 49, 0, 0),     Color(.25, .25, .75)), // pared der
        Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25)), // pared detras
        Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75)), // suelo
        Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25)), // techo
        // esferas normales para el resto de muestreos
        Sphere(16.5, Point(-23, -24.3, -34.6),  Color(.2, .3, .4)), 
        Sphere(16.5, Point(23, -24.3, -3.6),     Color(.4, .3, .2)),

        // para los valores del aluminio el 0.9 es para un gris claro que no afecta. para el oro es un amarillo para que de como oro

        // esferas de aluminio y oro para microfacet
       /* Sphere(16.5, Point(-23, -24.3, -34.6),  Color(.9, .9, .9), Color(), CONDUCTOR, 0.3, 
           Color(1.44, 0.96, 0.61), Color(7.47, 6.52, 5.29)), // Aluminio
        Sphere(16.5, Point(23, -24.3, -3.6),     Color(.9, .8, .3), Color(), CONDUCTOR, 0.3,
           Color(0.143, 0.374, 1.442), Color(3.982, 2.386, 1.603)), // Oro */

        // fuentes luminosas   
        //FuenteLuminosa(10.5, Point(0, 24.3, 0), Color(1, 1, 1)) // fuente normal
    	//FuentePuntual(Point(0, 24.3, 0)) // fuente puntual 

        FuenteLuminosa(10.5, Point(-23, 24.3, 0), Color(1, 1, 1), Color(12, 5, 5)),
        FuenteLuminosa(5, Point(23, 24.3, -50), Color(1, 1, 1), Color(5, 5, 12)) 
};

Sphere& fuenteLuminosa = spheres[N_ARREGLO]; // se define el arreglo de esferas.

// limita el valor de x a [0,1]
inline double clamp(const double x) { 
	if(x < 0.0)
		return 0.0;
	else if(x > 1.0)
		return 1.0;
	return x;
}

// convierte un valor de color en [0,1] a un entero en [0,255]
inline int toDisplayValue(const double x) {
	return int( pow( clamp(x), 1.0/2.2 ) * 255 + .5); 
}

// calcular la intersección del rayo r con todas las esferas
// regresar true si hubo una intersección, falso de otro modo
// almacenar en t la distancia sobre el rayo en que sucede la interseccion
// almacenar en id el indice de spheres[] de la esfera cuya interseccion es mas cercana
inline bool intersect(const Ray &r, double &t, int &id) {
	t = 1e20; //Un número muy grande 
	int numSpheres = N_ESFERAS; //El número de las esferas
	for(int i = 0; i < numSpheres; i++){ //verificamos todas las esferas
		double dist = spheres[i].intersect(r); //guardamos la distancia en la que interseccta la esfera i
		if(dist > 1e-4 && dist < t) // si la distancia es mayor que ~0 y menor que t 
		{
			t = dist; // ahora t es la menor encontrada
			id = i; // y se guarda la id
		}
	}
	return t < 1e19; // Si se encontro al menos una, regresa verdadero
}

// se usa en mapeos de texturas, cálculo de brdf, y para pasar de local a globall.
// todos los vectores son perpendiculares entre sí
void coordinateSystem(Vector &n, Vector &s, Vector &t){ // es una fucnión para definir el sistema de coordenadas en 179
	if (fabs(n.x) > fabs(n.y)) // el verdadero vector que nos dan es n, el resto los tenemos que calcular. si x es mayor a y
	{ // además, n debe de ya estar normalizado
        //par cálcular t, vemos el eje de n que domina más, para 
		float invLen = 1.0f / sqrt(n.x * n.x + n.z * n.z); // entonces sacamos un invLen el cual escala el vector para que t tenga longitud 1.
		t = Vector(n.z * invLen, 0.0f, -n.x * invLen); // y calulamos el vector t garantizando que t y n sean perpendiculares
	} else {
		float invLen = 1.0f / sqrt(n.y * n.y + n.z * n.z);
		t = Vector(0.0f, n.z * invLen, -n.y * invLen);
	}
	s = t%n; // en base al vector t y n, se calcula un vector ortogonal
}

Vector esfericasACartesianas(double cosTheta, double sinTheta, double phi){ // declaramos una función para pasar de eesféricas a cartecianas. las fóumulas están en el pdf del tema 2
	return Vector(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta); // podemos ignorar el radio ya que el vector resultante ya está normalizado
}
// --------------------------- HORA 2 ------------------------------------------------------------------------

// es una fución para cambiar de local a global en 180. Recibe el vector local y los otros vectores los cuales son los ejes globales
// para pasar de local a global se pasa por una matriz M, y se multiplica por el vector local
// esta fómula también es del tema 2
// se usan en el muestreo hemisferico, coseno y microfacet. 
Vector localAGlobal(const Vector &local, Vector &n, Vector &s, Vector &t){  
	return Vector( // regresa un vector el cual ya está en coordenadas globales
		s.x * local.x + t.x * local.y + n.x * local.z, // x es la suma de los componentes en x de los vectores globales por los x, y, z de los locales
		s.y * local.x + t.y * local.y + n.y * local.z, // y es la suma de los componentes y de los globales por los x, y, z de los localaes
		s.z * local.x + t.z * local.y + n.z * local.z  // y lo mismo para z
	);
}

// es una función para cambiar de  global a local en 189.  Recibe el vector global y los vectores del sistema local a cambiar
// aquí hay una M a la -1, por la cual vamos a multiplicar nuestro vector globlal
// se usa en la brdf de midrofacet, en el término g y g1
Vector globalALocal(const Vector &global, Vector &n, Vector &s, Vector &t){
	return Vector( // y regresa mos el vector pasado a las coordenadas locales
		global.dot(s), // se usa el producot punto ya que al hacer la multiplicación por 
		global.dot(t), // la matriz M a la -1 nos daba los componentes de x con x y y con y de la 
		global.dot(n)  // la suma de cada componente como con el producto punto
	);
}

// se hace el muestreo hemisférico en 197 que devuelve una estructura donde se almacenan los resultanos
// necesrios pasa montecarlo. Escogemos un punto aleatorio uniforme en el hemisferio y de allí vemos todas 
// las contribuciones que puede tener de todos los lados
// recibe un punto aleatorio del monteCarlo, un vector n igual dado por monteCarlo y la probabilidad que
// se calculará en esta función
// las fórmulas están en el tema 3
MuestreoResult MuestreoUniformeHemisferico(const Point &x, const Vector &n, double &prob) {
	double u1 = (double)rand() / RAND_MAX; // primero definimos a los números aleatorios
	double u2 = (double)rand() / RAND_MAX; // deben ser dos distintos

	double theta = acos(u1); // el ángulo polar para sólo muestrear un hemisferio
	double phi = 2.0 * M_PI * u2; // ángulo azimutal también sólo para un hemisferio
    //con estas fómulas ya tenemos elegidos nuestros calores para un punto aleatorio uniforme en el hemisferio
    // necesitamos estas fómulas y no sólo los números aleatorios ya que si lo hacemos así no sería uniforme
    // si vemos un globo terraqueo los cuadritos dearriba son más chikitos y los del ecuador más grandes.
    // Habría más probabilidad de caer en el acuador que en los polos.
    // entonces, con estas coordendas esféricas, las tenemso que pasar a cartesianas, las cuales están en 
    // local y de allí se pasarán a global
	double cosTheta = cos(theta);  // esto es para hacer el sistema de coordenadas con la normal de monteCarlo
	double sinTheta = sin(theta);  // igual 

	Vector s, t, normal = n; // delclaramos los vectores para hacer el sistema de coordenadas local
	coordinateSystem(normal, s, t); // y creamos el sistema de coordenadas tomando la normal que nos mandó monteCarlo 

	Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi); //las coordenadas que tenemos son esféricas y hacemos el vector de dirección
	Vector wi = localAGlobal(localDir, normal, s, t); // las pasamos a global, que sería nuestro vector wi 
	wi = wi.normalize(); //y lo normalizamos por que sólo nos interesa la dirección
	 
	prob = 1.0 / (2.0 * M_PI); // probabilidad uniforme sobre el hemisferio

	return MuestreoResult(wi, prob, Point()); // y regresamos la dirección y la probabilidad. Como este muestreo no va a una 
                                              // fuent de luz específica, entonces ponemos que es en 0 
}

// --------------------------- HORA 3 ------------------------------------------------------------------------


MuestreoResult MuestreoCosenoHemisferico(const Point &x, const Vector &n, double &prob) {
	double u1 = (double)rand() / RAND_MAX;
	double u2 = (double)rand() / RAND_MAX;

	double theta = asin(sqrt(u1));  
	double phi = 2.0 * M_PI * u2;

	double cosTheta = cos(theta); // para las coordenadas locales 
	double sinTheta = sin(theta);

	Vector s, t, normal = n;
	coordinateSystem(normal, s, t);

	Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi);
    Vector wi = localAGlobal(localDir, normal, s, t);
    wi = wi.normalize();

    prob = cosTheta / M_PI;
    
    return MuestreoResult(wi, prob, Point());
}

MuestreoResult MuestreoFuentePuntual(const Point &x, double &prob) {
	Vector wi = fuenteLuminosa.p - x; // apunta al centro de la fuente
	prob = 1.0;
	return MuestreoResult(wi, prob, fuenteLuminosa.p);
}
// --------------------------- HORA 4 ------------------------------------------------------------------------

Point puntoAleatorioEnEsfera(const Point &centro, double radio) {
	// vomo en muestreo hemisférico
	double u1 = (double)rand() / RAND_MAX;
	double u2 = (double)rand() / RAND_MAX;
	double theta = acos(2.0 * u2 - 1.0);
	double phi = 2.0 * M_PI * u1;

	double cosTheta = cos(theta);
    double sinTheta = sin(theta);

	Vector puntoUnitario = esfericasACartesianas(cosTheta, sinTheta, phi);
    return centro + puntoUnitario * radio;
}

MuestreoResult MuestreoArea(const Point &x, double &prob) {
    Point puntoEnFuente = puntoAleatorioEnEsfera(fuenteLuminosa.p, fuenteLuminosa.r);

    Vector wi = puntoEnFuente - x;
    double distanciaAlCuadrado = wi.dot(wi);
    wi = wi.normalize();

    double areaEsfera = 4.0 * M_PI * fuenteLuminosa.r * fuenteLuminosa.r;
    Vector normalFuente = (puntoEnFuente - fuenteLuminosa.p).normalize(); // normal en el punto de la fuente 
    Vector direccionHaciaX = (x - puntoEnFuente).normalize(); // dirección desde el punto en la fuente hacia el punto x
    double cosThetaFuente = normalFuente.dot(direccionHaciaX);  // coseno del ángulo entre la normal de la fuente y la dirección hacia x
     
    // si el punto está en la parte de atrás de la fuente la probabiildad es baja
    if (cosThetaFuente <= 0) {
        cosThetaFuente = 1e-6;
    }

    prob = distanciaAlCuadrado / (areaEsfera * fabs(cosThetaFuente)); // = 1/area * (r^2 / |cos(theta_fuente)|)
    return MuestreoResult(wi, prob, puntoEnFuente);
}

// --------------------------- HORA 5 ------------------------------------------------------------------------

// ***************** LUNES ******************************************************************************
MuestreoResult MuestreoAnguloSolido(const Point &x, double &prob) {
    Vector dirCentro = fuenteLuminosa.p - x; 
    double distancia = sqrt(dirCentro.dot(dirCentro)); // distancia del punto x al centro de la fuente
   
    double sinThetaMax = fuenteLuminosa.r / distancia; 
    double thetaMax = asin(sinThetaMax);
    
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    
    double cosTheta = (1.0 - u1) + u1 * cos(thetaMax);
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    double phi = 2.0 * M_PI * u2;

    Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi);

    Vector centroNormalizado = dirCentro.normalize();
    Vector s, t;
    coordinateSystem(centroNormalizado, s, t);

    Vector wi = localAGlobal(localDir, centroNormalizado, s, t);
    wi = wi.normalize();

    //punto de intersección en la fuente
    Ray rayo(x, wi);
    double tInterseccion = fuenteLuminosa.intersect(rayo);
    Point puntoEnFuente = x + wi * tInterseccion;

    prob = 1.0 / (2.0 * M_PI * (1.0 - cos(thetaMax)));
    
    return MuestreoResult(wi, prob, puntoEnFuente);
}
// --------------------------- HORA 6 ------------------------------------------------------------------------

double xMas(double x) {
    return x > 0 ? 1.0 : 0.0;
}

double DistribucionBeckmann(const Vector &wh, double aspereza) {
    double cosTheta = wh.z; // ángulo con la normal
    double alpha = aspereza;
    double cos2Theta = cosTheta * cosTheta;
    double tan2Theta = (1.0 - cos2Theta) / cos2Theta;
    
    double chi = xMas(cosTheta);
    return chi * exp(-tan2Theta / (alpha * alpha)) / 
           (M_PI * alpha * alpha * cos2Theta * cos2Theta);
}

double TerminoG1(const Vector &w, const Vector &wh, double aspereza) {
    double cosTheta = w.z;
    if (cosTheta <= 0) return 0.0;
    
    double alpha = aspereza;
    double cos2Theta = cosTheta * cosTheta;
    double tan2Theta = (1.0 - cos2Theta) / cos2Theta;
    double tanTheta = sqrt(tan2Theta);
    
    double a = 1.0 / (alpha * tanTheta);
    if (a >= 1.6) return 1.0;
    
    double a2 = a * a;
    return (3.535 * a + 2.181 * a2) / (1.0 + 2.276 * a + 2.577 * a2);
}
// --------------------------- HORA 7 ------------------------------------------------------------------------

double TerminoG(const Vector &wo, const Vector &wi, const Vector &wh, double aspereza) {
    return TerminoG1(wo, wh, aspereza) * TerminoG1(wi, wh, aspereza);
}

Color FresnelConductor(double cosTheta, const Color &eta, const Color &k) {
    double cos2Theta = cosTheta * cosTheta;
    double sin2Theta = 1.0 - cos2Theta;
    
    Color eta2 = Color(eta.x * eta.x, eta.y * eta.y, eta.z * eta.z);
    Color k2 = Color(k.x * k.x, k.y * k.y, k.z * k.z);
    
    Color t0 = eta2 - k2 - Color(sin2Theta, sin2Theta, sin2Theta);
    Color a2Masb2 = Color(sqrt(t0.x * t0.x + 4.0 * eta2.x * k2.x),
                          sqrt(t0.y * t0.y + 4.0 * eta2.y * k2.y),
                          sqrt(t0.z * t0.z + 4.0 * eta2.z * k2.z));
    
    Color t1 = a2Masb2 + Color(cos2Theta, cos2Theta, cos2Theta);
    Color a = Color(sqrt(0.5 * (a2Masb2.x + t0.x)),
                   sqrt(0.5 * (a2Masb2.y + t0.y)),
                   sqrt(0.5 * (a2Masb2.z + t0.z)));
    
    Color t2 = Color(2.0 * cosTheta * a.x, 2.0 * cosTheta * a.y, 2.0 * cosTheta * a.z);
    Color Rs = Color((t1.x - t2.x) / (t1.x + t2.x),
                    (t1.y - t2.y) / (t1.y + t2.y),
                    (t1.z - t2.z) / (t1.z + t2.z));
    
    Color t3 = Color(cos2Theta * a2Masb2.x + sin2Theta * sin2Theta,
                    cos2Theta * a2Masb2.y + sin2Theta * sin2Theta,
                    cos2Theta * a2Masb2.z + sin2Theta * sin2Theta);
    
    Color t4 = t2 * sin2Theta;
    Color Rp = Color(Rs.x * (t3.x - t4.x) / (t3.x + t4.x),
                    Rs.y * (t3.y - t4.y) / (t3.y + t4.y),
                    Rs.z * (t3.z - t4.z) / (t3.z + t4.z));
    
    return Color(0.5 * (Rs.x + Rp.x), 0.5 * (Rs.y + Rp.y), 0.5 * (Rs.z + Rp.z));
}
// --------------------------- HORA 8 ------------------------------------------------------------------------

Vector MuestreoMicrofacet(double aspereza, double u1, double u2) {
    double alpha = aspereza;
    double tan2Theta = -alpha * alpha * log(1.0 - u1);  //  D(wh) = χ⁺(cos θ) * exp(-tan^2θ/α²) / (π α^2 cos^2θ)
    double cos2Theta = 1.0 / (1.0 + tan2Theta); // 1 + tan^2θ = sec^2θ = 1/cos^2θ
    double cosTheta = sqrt(cos2Theta); 

    double sin2Theta = 1.0 - cos2Theta; // sin^2θ + cos^2θ = 1  =>  sin^2θ = 1 - cos^2θ
    double sinTheta = sqrt(sin2Theta); 
    double phi = 2.0 * M_PI * u2;

    return esfericasACartesianas(cosTheta, sinTheta, phi);
}

Color BRDFMicrofacet(const Sphere &obj, const Vector &wi, const Vector &wo, const Vector &n) {
    Vector s, t, normal = n;
    coordinateSystem(normal, s, t);
    
    Vector wiLocal = globalALocal(wi, normal, s, t);
    Vector woLocal = globalALocal(wo, normal, s, t);
    
    if (wiLocal.z <= 1e-6 || woLocal.z <= 1e-6) return Color();
    
    Vector whLocal = (wiLocal + woLocal).normalize();
    
    double D = DistribucionBeckmann(whLocal, obj.aspereza);
    double G = TerminoG(woLocal, wiLocal, whLocal, obj.aspereza);
    Color F = FresnelConductor(abs(wiLocal.dot(whLocal)), obj.eta, obj.kappa); 
    
    double denominador = 4.0 * abs(woLocal.z) * abs(wiLocal.z); 
    if (denominador < 1e-6) return Color();
    
    return Color(F.x * D * G / denominador, F.y * D * G / denominador, F.z * D * G / denominador);
}
// --------------------------- HORA 8 ------------------------------------------------------------------------

// función general de Monte Carlo que acepta cualquier estrategia de muestreo
Color monteCarlo(int N, Vector &n, Point &x, const Sphere &obj, 
                 MuestreoResult (*estrategiaMuestreo)(const Point&, const Vector&, double&)) {
    Color sum;
    Point puntoSeguro = x + n * 1e-4;

    for (int i = 0; i < N; ++i) {
        double prob;
        MuestreoResult resultado = estrategiaMuestreo(puntoSeguro, n, prob);
        Vector wi = resultado.wi;

        Ray shadowRay(puntoSeguro, wi);
        double t_shadow;
        int id_shadow;

        bool hayInterseccion = intersect(shadowRay, t_shadow, id_shadow);

        if (!hayInterseccion) {
            continue;
        }
        
        const Sphere &objIntersectado = spheres[id_shadow];
        Color Le = objIntersectado.e;
        Color fr = obj.c * (1.0 / M_PI);
        double cosTheta = n.dot(wi);

        if (cosTheta > 0) {
            Color contribucion = Le.mult(fr) * (cosTheta / prob);
            sum = sum + contribucion; 
        }
    }
    return sum * (1.0 / N);
}
// --------------------------- HORA 9 ------------------------------------------------------------------------

Color monteCarloHemisfer(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoUniformeHemisferico);
}

Color monteCarloCoseno(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoCosenoHemisferico);
}

MuestreoResult MuestreoAreaWrapper(const Point &x, const Vector &n, double &prob) {
    return MuestreoArea(x, prob);
}

Color monteCarloArea(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoAreaWrapper);
}

MuestreoResult MuestreoAnguloSolidoWrapper(const Point &x, const Vector &n, double &prob) {
    return MuestreoAnguloSolido(x, prob);
}

Color monteCarloAnguloSolido(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoAnguloSolidoWrapper);
}
// --------------------------- HORA 10 ------------------------------------------------------------------------

Color monteCarloFuentePuntual(Vector &n, Point &x, const Sphere &obj) {
    Point puntoSeguro = x + n * 1e-4;
    
    // obtener dirección y punto
    double prob;
    MuestreoResult resultado = MuestreoFuentePuntual(puntoSeguro, prob);
    Vector wi = resultado.wi;
    Point puntoEnFuente = resultado.puntoEnFuente;
    
    // calcular distancia a la fuente
    Vector direccionCompleta = wi;
    double distanciaAlCuadrado = direccionCompleta.dot(direccionCompleta);
    double distancia = sqrt(distanciaAlCuadrado);
    wi = wi.normalize();
    // verificar si hay obstáculos 
    Ray shadowRay(puntoSeguro, wi);
    double t_shadow;
    int id_shadow;
    
    if (intersect(shadowRay, t_shadow, id_shadow)) {
        // si intersecta algo antes de llegar a la fuente, está en sombra
        if (t_shadow < distancia - 1e-4) {
            return Color(); // en sombra
        }
    }

    Color Le = fuenteLuminosa.e;
    Color fr = obj.c * (1.0 / M_PI); // BRDF lambertiana
    double cosTheta = n.dot(wi);
    
    if (cosTheta > 0) {
        // para fuente puntual, la intensidad decrece con 1/r^2
        Color contribucion = Le.mult(fr) * (cosTheta / distanciaAlCuadrado);
        return contribucion;
    }
    
    return Color(); // no hay contribución si cosTheta <= 0
}
// --------------------------- HORA 11 ------------------------------------------------------------------------

Color pathTracing(const Ray &ray, int profundidad, int maxProfundidad) {
    if (profundidad >= maxProfundidad) return Color();
    
    double t;
    int id;
    if (!intersect(ray, t, id)) return Color();
    
    const Sphere &obj = spheres[id];
    Point x = ray.o + ray.d * t;
    Vector n = (x - obj.p).normalize();
    
    // Si es fuente luminosa, retornar emisión
    if (obj.e.x > 0 || obj.e.y > 0 || obj.e.z > 0) {
        return obj.e;
    }
 
    Point puntoSeguro = x + n * 1e-4;
    Color L(0,0,0);
    
    // Iluminación directa. Se puede cambiar el tipo de muestreo
    L = L + monteCarloCoseno(N_,n, x, obj);
    
    // Muestreo según material
    Vector wo = Vector(0,0,0) - ray.d;
    
    if (obj.material == DIFFUSO) {
        double prob;
        MuestreoResult resultado = MuestreoCosenoHemisferico(puntoSeguro, n, prob);
        Vector wi = resultado.wi;
// --------------------------- HORA 12 ------------------------------------------------------------------------
        
        if (prob > 1e-6 && n.dot(wi) > 1e-6) {
            Ray nuevoRayo(puntoSeguro, wi);
            Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad);
            Color fr = obj.c * (1.0 / M_PI);
            double cosTheta = n.dot(wi);
            
            double factor = (cosTheta / prob);
            L = L + Color(fr.x * Li.x * factor, fr.y * Li.y * factor, fr.z * Li.z * factor);
        }
    } 
    else if (obj.material == CONDUCTOR) {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = (double)rand() / RAND_MAX;
        
        Vector s, t;
        coordinateSystem(n, s, t);
        
        Vector whLocal = MuestreoMicrofacet(obj.aspereza, u1, u2);
        Vector wh = localAGlobal(whLocal, n, s, t);
        wh = wh.normalize();
        
        Vector wi = wh * (2.0 * wo.dot(wh)) - wo;
        wi = wi.normalize();
        
        if (n.dot(wi) > 1e-6) {
            double D = DistribucionBeckmann(whLocal, obj.aspereza);
            double G = TerminoG(globalALocal(wo, n, s, t), globalALocal(wi, n, s, t), whLocal, obj.aspereza);
            Color F = FresnelConductor(fabs(wo.dot(wh)), obj.eta, obj.kappa);
            
            double probMicrofacet = D * whLocal.z / (4.0 * fabs(wo.dot(wh)));
            double denom = probMicrofacet * fabs(n.dot(wh));
// --------------------------- HORA 12 ------------------------------------------------------------------------

// ************************************* MARTES **************************************************
            if (denom > 1e-6) {
                Ray nuevoRayo(puntoSeguro, wi);
                Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad);
                double factor = (G * fabs(wo.dot(wh)) / denom);
                
                L = L + Color(F.x * Li.x * factor, F.y * Li.y * factor, F.z * Li.z * factor);
            }
        }
    }
    
    return Color(fmin(fmax(L.x, 0.0), 1.0), fmin(fmax(L.y, 0.0), 1.0), fmin(fmax(L.z, 0.0), 1.0));
}


bool shade(const Ray &r, Point &x, Vector &n, const Sphere *&obj) {
    double t;
    int id = 0;
    
    // determinar que esfera (id) y a que distancia (t) el rayo intersecta
    if (!intersect(r, t, id))
        return false;  // el rayo no intersecto objeto
    
    obj = &spheres[id];
    x = r.o + r.d * t; // punto de intersección
    n = (x - obj->p).normalize(); // normal en el punto de intersección
    
    return true;
}
// --------------------------- HORA 13 ------------------------------------------------------------------------

Color shadeHemisfer(const Ray &r) {
    Point x;
    Vector n;
    const Sphere *obj;
    
    if (!shade(r, x, n, obj))
        return Color();  // el rayo no intersecto objeto, return Vector() == negro
    // Si el objeto emite luz, retornar directamente la emisión
    if (obj->e.x > 0 || obj->e.y > 0 || obj->e.z > 0) { 
        return obj->e; 
    }
    return obj->e + monteCarloHemisfer(N_, n, x, *obj);
}

Color shadeCoseno(const Ray &r) {
    Point x;
    Vector n;
    const Sphere *obj;
    
    if (!shade(r, x, n, obj))
        return Color();  // el rayo no intersecto objeto, return Vector() == negro
    // Si el objeto emite luz, retornar directamente la emisión
    if (obj->e.x > 0 || obj->e.y > 0 || obj->e.z > 0) { 
        return obj->e; 
    }
    return obj->e + monteCarloCoseno(N_, n, x, *obj);
}

Color shadeFuentePuntual(const Ray &r) {
    Point x;
    Vector n;
    const Sphere *obj;
// --------------------------- HORA 13 ------------------------------------------------------------------------
    
    if (!shade(r, x, n, obj))
        return Color();  // el rayo no intersecto objeto, return Vector() == negro
    // Si el objeto emite luz, retornar directamente la emisión
    if (obj->e.x > 0 || obj->e.y > 0 || obj->e.z > 0) { 
        return obj->e; 
    }
    return obj->e + monteCarloFuentePuntual(n, x, *obj);
}

Color shadeArea(const Ray &r) {
    Point x;
    Vector n;
    const Sphere *obj;
    
    if (!shade(r, x, n, obj))
        return Color();  // el rayo no intersecto objeto, return Vector() == negro
    // Si el objeto emite luz, retornar directamente la emisión
    if (obj->e.x > 0 || obj->e.y > 0 || obj->e.z > 0) { 
        return obj->e; 
    }
    return obj->e + monteCarloArea(N_, n, x, *obj);
}

Color shadeAnguloSolido(const Ray &r) {
    Point x;
    Vector n;
    const Sphere *obj;
    
    if (!shade(r, x, n, obj))
        return Color();  // el rayo no intersecto objeto, return Vector() == negro
    // si el objeto emite luz, retornar directamente la emisión
    if (obj->e.x > 0 || obj->e.y > 0 || obj->e.z > 0) { 
        return obj->e; 
    }
    return obj->e + monteCarloAnguloSolido(N_, n, x, *obj);
}
// --------------------------- HORA 14 ------------------------------------------------------------------------

Color shadeMicrofacet(const Ray &r) {
    Color L(0,0,0);
    int spp = 64; // samples per pixel
    for (int i = 0; i < spp; ++i) {
        Color sample = pathTracing(r, 0, 8);
        // acumulación
        L = Color(L.x + sample.x / spp,
                 L.y + sample.y / spp,
                 L.z + sample.z / spp);
    }
    return L;
}

int main(int argc, char *argv[]) {
	int w = 1024, h = 768; // image resolution
  
	// fija la posicion de la camara y la dirección en que mira
	Ray camera( Point(0, 11.2, 214), Vector(0, -0.042612, -1).normalize() );
	// parametros de la camara
	Vector cx = Vector( w * 0.5095 / h, 0., 0.); 
	Vector cy = (cx % camera.d).normalize() * 0.5095;
  
	// auxiliar para valor de pixel y matriz para almacenar la imagen
	Color *pixelColors = new Color[w * h];

	srand(time(NULL));

	#pragma omp parallel for
	for(int y = 0; y < h; y++) 
	{ 
// --------------------------- HORA 15 ------------------------------------------------------------------------

		// recorre todos los pixeles de la imagen
		fprintf(stderr,"\r%5.2f%%",100.*y/(h-1));
		for(int x = 0; x < w; x++ ) {
			int idx = (h - y - 1) * w + x; // index en 1D para una imagen 2D x,y son invertidos
			Color pixelValue = Color(); // pixelValue en negro por ahora
			// para el pixel actual, computar la dirección que un rayo debe tener
			Vector cameraRayDir = cx * ( double(x)/w - .5) + cy * ( double(y)/h - .5) + camera.d;
			
			// computar el color del pixel para el punto que intersectó el rayo desde la camara

			//Descomentar el tipo de muestreo que se quiere usar. Recordar también descomentar la fuente que se va a usar, normal p puntual
			//pixelValue = shadeHemisfer( Ray(camera.o, cameraRayDir.normalize()) ); // Hemisferico
			//pixelValue = shadeCoseno( Ray(camera.o, cameraRayDir.normalize()) ); // Coseno hemisferico
			//pixelValue = shadeFuentePuntual( Ray(camera.o, cameraRayDir.normalize()) ); // Fuente Puntual
			//pixelValue = shadeArea( Ray(camera.o, cameraRayDir.normalize()) ); // Area
		    //pixelValue = shadeAnguloSolido( Ray(camera.o, cameraRayDir.normalize()) ); // Ángulo Sólidp
			pixelValue = shadeMicrofacet( Ray(camera.o, cameraRayDir.normalize()) ); // Microfacet
			

			// limitar los tres valores de color del pixel a [0,1]
			pixelColors[idx] = Color(clamp(pixelValue.x), clamp(pixelValue.y), clamp(pixelValue.z));
		}
	}

	fprintf(stderr,"\n");
// --------------------------- HORA 16 ------------------------------------------------------------------------

	FILE *f = fopen("image.ppm", "w");
	// escribe cabecera del archivo ppm, ancho, alto y valor maximo de color
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int p = 0; p < w * h; p++) 
	{ // escribe todos los valores de los pixeles
    		fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), 
				toDisplayValue(pixelColors[p].z));
  	}
  	fclose(f);

  	delete[] pixelColors;

	return 0;
}
// --------------------------- HORA 17 ------------------------------------------------------------------------
// ****************************** Martes en la madrugada repasar preguntas () ******************************************