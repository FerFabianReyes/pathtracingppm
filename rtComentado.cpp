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

// El muestreo de conseno hemisférico es una mejora al muestreo de hemisferio.
// recordando que en el muestreo del heisferio le damos la misma importancia a todos los puntos del hemisferio
// ahora le damos más importancia a la normal que "sale" de la esfera, ya que entre más 
// recto esté con la luz, más contribución tendra, en lugar de las que están inclinadas a la fuente de luz
// estas fórmulas están en el tema 3
MuestreoResult MuestreoCosenoHemisferico(const Point &x, const Vector &n, double &prob) {
	double u1 = (double)rand() / RAND_MAX; // igual necesitamos dos números aleatorios
	double u2 = (double)rand() / RAND_MAX;

	double theta = asin(sqrt(u1)); // cálculamos el ángulo polar
	double phi = 2.0 * M_PI * u2; // y el azimutal
    // con estas fórmulas, tenemos nuestras coordenadas en un punto que favorece a lo descrito arriba
    // entre más cerca de la normal esté, mayor contribución, entre más inclinado, menor contribución 

	double cosTheta = cos(theta); // igual que en el muestreo de hemisferio, necesitamos pasar a  
	double sinTheta = sin(theta); // otro sistema de coordenadas.

	Vector s, t, normal = n;
	coordinateSystem(normal, s, t); // creamos el sistema de coordenadas local

	Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi); // pasamos nuestras coordenadas cartesianas a locales
    Vector wi = localAGlobal(localDir, normal, s, t); // de locales a globales
    wi = wi.normalize(); // y normalizamos ya que sólo nos interesa la dirección

    prob = cosTheta / M_PI; // ahora nuestra probabilidad ya no es constante, depende del ángulo polar
    //  se divide entre π porque es la constante de normalización para que la integral de la PDF sobre el hemisferio sea igual a 1.
    
    return MuestreoResult(wi, prob, Point()); // y regresamos la dirección y la probabilidad. Como este muestreo no va a una 
                                              // fuent de luz específica, entonces ponemos que es en 0 
}


// Para el muestreo de una fuente puntual, debemos de conocer su posición, ya que como sólo es un punto, atinarle es
// imposible. En este muestreo conocemos la posición de la fuente de luz. 
MuestreoResult MuestreoFuentePuntual(const Point &x, double &prob) { // recibimos un punto y la probabilidad a calcular
	Vector wi = fuenteLuminosa.p - x; // apunta al centro de la fuente
	prob = 1.0; // la probabilidad es 1 por que siempre le daremos a la fuente
	return MuestreoResult(wi, prob, fuenteLuminosa.p);  // y damos el resultado. Ahora como punto ponemos el punto de la fuente de luz
    // no normalizamos el vector wi por que en el monteCarlo de la fuete puntual necesitamos sacara la distancia
}
// --------------------------- HORA 4 ------------------------------------------------------------------------

// esta función sólo nos escoge un punto aleatorio en la superficie de la fuente de luz en 248.
// se usa nada más para el muestreo de área
// recibe el centro de la fuente y el radio
Point puntoAleatorioEnEsfera(const Point &centro, double radio) {
	// en el tema 4 se dice que se usa como en muestreo hemisférico
	double u1 = (double)rand() / RAND_MAX; // necesitamos dos números aleatorios
	double u2 = (double)rand() / RAND_MAX;
	double theta = acos(2.0 * u2 - 1.0); // con 2.0 * u2 - 1.0 nos aseguramos que sólo esté entre -1 y 1 el ángulo polar
	double phi = 2.0 * M_PI * u1;  // ángulo azimutal

	double cosTheta = cos(theta); // igual para el sistema de coordenadas local
    double sinTheta = sin(theta); 

	Vector puntoUnitario = esfericasACartesianas(cosTheta, sinTheta, phi); // las coordenadas esféricas que tenemos las pasamos a cartesianas
    return centro + puntoUnitario * radio; // ahora a este punto lo multiplicamos para que esté en el borde de la 
                                           // fuente de luz y lo trazladamos al centro 
}


// para el muestreo de área en 262 es que tenemos una fuente de luz que sabemos donde está
// es como el muestreo hemisférico pero en lugar de buscar la fuente de luz, hacemos
// el muestreo en el área de la fuente de luz
// las formulas están en el tema 4
MuestreoResult MuestreoArea(const Point &x, double &prob) { // recibimos un punto y la variable a calcular de la probabilidad
    Point puntoEnFuente = puntoAleatorioEnEsfera(fuenteLuminosa.p, fuenteLuminosa.r); // escogemos un punto aleatorio en la fuente

    Vector wi = puntoEnFuente - x; // cáclulamos el vector desde el punto x hacia la fuente luminosa
    double distanciaAlCuadrado = wi.dot(wi);
    wi = wi.normalize(); // y normalizamos por que sólo nos interesa la dirección

    double areaEsfera = 4.0 * M_PI * fuenteLuminosa.r * fuenteLuminosa.r; // el áre a de 4*pi*r^2
    Vector normalFuente = (puntoEnFuente - fuenteLuminosa.p).normalize(); // normal en el punto de la fuente 
    Vector direccionHaciaX = (x - puntoEnFuente).normalize(); // dirección desde el punto en la fuente hacia el punto x
    double cosThetaFuente = normalFuente.dot(direccionHaciaX);  // coseno del ángulo entre la normal de la fuente y la dirección hacia x
     
    // si el punto está en la parte de atrás de la fuente la probabiildad es baja
    // Es el coseno del ángulo entre la normal de la fuente y la dirección hacia el punto
    // Representa cuánta luz "escapa" de la fuente en esa dirección
    // Si cosThetaFuente = 0: la fuente está de lado (poca contribución)
    // Si cosThetaFuente < 0: la fuente está de espaldas (no hay contribución)
    if (cosThetaFuente <= 0) {
        cosThetaFuente = 1e-6;
    }

    prob = distanciaAlCuadrado / (areaEsfera * fabs(cosThetaFuente)); // = 1/area * (r^2 / |cos(theta_fuente)|)
    // la fórmula de la probabilidad viene de que de área se pasa a ángulo sólido. 
    // Es como convertir "probabilidad por metro cuadrado de la fuente" a "probabilidad por dirección que puedo mirar"
    //dA = área infinitesimal en la fuente
    //r² = la "dilución" de la luz con la distancia
    //cos(θ) = proyección (solo la componente perpendicular contribuye)
    //dΩ = ángulo sólido visto desde el punto

    return MuestreoResult(wi, prob, puntoEnFuente);
}

// --------------------------- HORA 5 ------------------------------------------------------------------------

// ***************** LUNES ******************************************************************************

// El muestreo de ángulo sólido en 283 es una mejora en el muestreo de área, ya que ahora se muestreoa directamente con
// un ángulo solido. Ahora, las muestras se obtienen en un cono que va desde el punto x hasta la fuente
// dentro de ese cono se muestrea
MuestreoResult MuestreoAnguloSolido(const Point &x, double &prob) { // recibimos un punto y probabilidad a calcular
    Vector dirCentro = fuenteLuminosa.p - x; // tenemos un vector desde el punto X hacia el centro de la fuente luminosa
    double distancia = sqrt(dirCentro.dot(dirCentro)); // distancia del punto x al centro de la fuente
   
    double sinThetaMax = fuenteLuminosa.r / distancia; 
    double thetaMax = asin(sinThetaMax); // thetaMax es el ángulo que delimita el borde del cono en el que estamos muestreando
    
    double u1 = (double)rand() / RAND_MAX; // necesitamos dos números aleatorios
    double u2 = (double)rand() / RAND_MAX; // para el muestreo
    
    double cosTheta = (1.0 - u1) + u1 * cos(thetaMax); // esto genera el rango que va desde cos(0) = 1 (dirección directa al centro) hasta cos(thetaMax) (borde del cono)
    double sinTheta = sqrt(1.0 - cosTheta * cosTheta);  // el seno se saca por identidades trigonométricas
    double phi = 2.0 * M_PI * u2;

    Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi); // las direcciones esféricas que obtuvimos las pasamos a cartesianas

    Vector centroNormalizado = dirCentro.normalize(); // normalizamos la dirección hacia el centro de la fuente para construir el sistema de coordenadas local
    Vector s, t; 
    coordinateSystem(centroNormalizado, s, t); // generamos el sistema de coordenadas con respecto a la fuente luminosa

    Vector wi = localAGlobal(localDir, centroNormalizado, s, t); // ahora lo pasamos a global
    wi = wi.normalize(); // y lo normalizamos

    //punto de intersección en la fuente
    Ray rayo(x, wi); // ahora creamos un rayo en la dirección muestreada
    double tInterseccion = fuenteLuminosa.intersect(rayo); // encontramos donde intersecta la fuente luminosa
    Point puntoEnFuente = x + wi * tInterseccion; //// se calcula el punto de la intersección

    prob = 1.0 / (2.0 * M_PI * (1.0 - cos(thetaMax))); // la probablididad es constante para todas las muestras dentro del cono
    
    return MuestreoResult(wi, prob, puntoEnFuente);
}
// --------------------------- HORA 6 ------------------------------------------------------------------------

// es una función en 316 que si x es mayor a 0 es 1 pero si no es 0
double xMas(double x) {
    return x > 0 ? 1.0 : 0.0;
}

// Es el término D para el microfacet 
// D representa la microsuperficie. Dice cuántas microfacetas están orientadas a wh
// Depende de eta que fija la esperaza del material
// Si D(wh) es alto en una dirección wh, significa que muchas microfacetas están alineadas con esa dirección
// Superficies rugosas (eta alto): D se dispersa, creando reflejos difusos y borrosos.
// Superficies lisas (eta bajo): D se concentra cerca de n, produciendo reflejos nítidos.
// Las ecuaciones están en el tema 5
double DistribucionBeckmann(const Vector &wh, double aspereza) {
    double cosTheta = wh.z; // z apunta a la normal de la superficie. wh es el vector de la microfaceta.
                        // En coordenadas locales, la normal es Vector(0,0,1)
                        // cos(θ) = wh · normal = (x,y,z) · (0,0,1) = z
    double alpha = aspereza; // definimos el valor de la aspereza
    double cos2Theta = cosTheta * cosTheta; // definimos cosTheta^2 pq después tiene que ser ^4. en las coordenadas esféricas wh.z = cos(θ)
    double tan2Theta = (1.0 - cos2Theta) / cos2Theta; // en la ecuación también nos piden que sea tan^2
                    //  tan²(θ) = sin²(θ) / cos²(θ).      sin²(θ) = 1 - cos²(θ)
    double x = xMas(cosTheta); // es una función de visibilidad o función de  máscara. 
                // χ = 1 si cosTheta > 0 (la microfaceta es visible)
                // χ = 0 si cosTheta ≤ 0 (la microfaceta no es visible)
                // es el ángulo entre el vector de dirección media y la normal
    return x * exp(-tan2Theta / (alpha * alpha)) / 
           (M_PI * alpha * alpha * cos2Theta * cos2Theta); //fórmuma para calcúlar el valor de B
}

// en la 331  calcula la fracción de luz que no es obstruida por las microfacetas en una dirección específica
double TerminoG1(const Vector &w, const Vector &wh, double aspereza) { // recibe un vector w, puede ser wi o wo, el vector wh y la asperza
    double cosTheta = w.z;  // es equivalente a w * n en coordenadas locales, ya que cos(θ) = wh · normal = (x,y,z) · (0,0,1) = z
    if (cosTheta <= 0) return 0.0; // si la contribución está derás de la superficie, no hay contribución
    // theta es el ángulo entre la dirección de la luz (wi o wo) y la normal (n).
    double alpha = aspereza; // definimos la aspereza
    double cos2Theta = cosTheta * cosTheta; // cos^2 cálcular tan
    double tan2Theta = (1.0 - cos2Theta) / cos2Theta; // tan(θ) = sin(θ)/cos(θ)
    double tanTheta = sqrt(tan2Theta); // y cálculamos tan para el término a. es la tangente del ángulo con la normal
    
    double a = 1.0 / (alpha * tanTheta);  // es el parámetro de rugosidad ajustado por ángulo 
    if (a >= 1.6) return 1.0; // 1.6 es un valor fijo definido en las fórmulas y si regresa 1 es por que entonces 
                              // no hay contribución, ya que *1 queda igual todo. 1.6 es un umbral empírico definido 
                              // por el modelo de Smith. Si regresa 1 es que no hay obstáculos de las microfacetas.
    
    double a2 = a * a; // si a fue menor que 1.6 entonces se regresa la otra fórmula. 
    return (3.535 * a + 2.181 * a2) / (1.0 + 2.276 * a + 2.577 * a2);  // fórmula del tema 5 para g1 para beckmann
}
// --------------------------- HORA 7 ------------------------------------------------------------------------

// se usa para ver cuánta luz podemos ver en realidad
// El término G Multiplica las obstrucciones para ambas direcciones (luz entrante y saliente).
// Ejemplo: Si G1(wi) = 0.5 y G1(wo) = 0.5 → G = 0.25 (solo el 25% de la luz se refleja).
// estamos usando coordenadas locales por que nos estamos poniendo en el lugar de la microfaceta
double TerminoG(const Vector &wo, const Vector &wi, const Vector &wh, double aspereza) { 
    return TerminoG1(wo, wh, aspereza) * TerminoG1(wi, wh, aspereza); // cálculamos la cantidad de luz que no esobstrida. por la entrada y la salida
}

// Es para la cantidad de energía que es reflejada en una superficie de la luz que entra
// reflectancía dependiente del ángulo. Fresnell ayuda con los colores de los materiales
// intensifica en los bordes, por eso es para materiales conductores
Color FresnelConductor(double cosTheta, const Color &eta, const Color &k) { 
    double cos2Theta = cosTheta * cosTheta;   
    double sin2Theta = 1.0 - cos2Theta;  // sigue siendo la identidad sin²(θ) = 1 - cos²(θ)
    // eta es el índice de refracción y kapa el coeficiente de extinción y thetha es el ángulo de incidencia
    Color eta2 = Color(eta.x * eta.x, eta.y * eta.y, eta.z * eta.z); // para hacerla al cuadrado sólo multiplicamos sus elementos
    Color k2 = Color(k.x * k.x, k.y * k.y, k.z * k.z); // igual para kappa
    
    Color t0 = eta2 - k2 - Color(sin2Theta, sin2Theta, sin2Theta); // sin²θ es un color ya que las operaciones se realizan en 
                                                                   // vectores y para poder restarlas ya que kappa y eta son colores
                                                                   // Pero la razón principal es por que al ponerlo como color podemos 
                                                                   // realizar las operaciones para los tres colores RGB de manera paralela
    Color a2Masb2 = Color(sqrt(t0.x * t0.x + 4.0 * eta2.x * k2.x), // esta fórmula igual está definida.sólo que está puesta en partes
                          sqrt(t0.y * t0.y + 4.0 * eta2.y * k2.y), // t0 =  η² - κ² - sin²θ y todo junto es  a² + b² = √[(η² - κ² - sin²θ)² + 4η²κ²]
                          sqrt(t0.z * t0.z + 4.0 * eta2.z * k2.z)); // igual se hace en cada parte del color para hacer los cálculos de manera simultanea
                                                                    // para los tres colores
    Color t1 = a2Masb2 + Color(cos2Theta, cos2Theta, cos2Theta); // t1 = a² + b² + cos²θ
    Color a = Color(sqrt(0.5 * (a2Masb2.x + t0.x)), // a = √[0.5*(a² + b² + η² - κ² - sin²θ)]
                   sqrt(0.5 * (a2Masb2.y + t0.y)), // igual para todas las líneas
                   sqrt(0.5 * (a2Masb2.z + t0.z))); // para el cálculo en cada color
    
    Color t2 = Color(2.0 * cosTheta * a.x, 2.0 * cosTheta * a.y, 2.0 * cosTheta * a.z); // t2 = 2a cosθ
    Color Rs = Color((t1.x - t2.x) / (t1.x + t2.x), // Rpe es reflectancia para la polarización perpendicular
                    (t1.y - t2.y) / (t1.y + t2.y), // se hace para todos los RGB
                    (t1.z - t2.z) / (t1.z + t2.z)); 
    
    Color t3 = Color(cos2Theta * a2Masb2.x + sin2Theta * sin2Theta, // (a² + b²)cos²θ + sin⁴θ
                    cos2Theta * a2Masb2.y + sin2Theta * sin2Theta, //igual para todo
                    cos2Theta * a2Masb2.z + sin2Theta * sin2Theta);
    
    Color t4 = t2 * sin2Theta; // 2a cosθ sin²θ
    Color Rp = Color(Rs.x * (t3.x - t4.x) / (t3.x + t4.x), // Ra es reflectancia para la polarización paralela
                    Rs.y * (t3.y - t4.y) / (t3.y + t4.y), // igual para todos
                    Rs.z * (t3.z - t4.z) / (t3.z + t4.z));
    
    return Color(0.5 * (Rs.x + Rp.x), 0.5 * (Rs.y + Rp.y), 0.5 * (Rs.z + Rp.z)); // y regresamos el promedio de  las dos para cada color
    // Rpe y Rpa dependen de las propiedades eta y kappa. si eta es alto, aumenta la reflectancia, si kappa es alto absorve más
    // luz pero tambien refleja más. Ambas son la parte real e imaginaria del índice de refracción complejo. 
    // eta, la parte real, controla la velocidad de fase de la luz. A mayor eta, aumenta la reflectancia y más en los bordes
    // kappa es la parte imaginaria, cuantifica la aoborción de la luz en el material. Un kapa alto indica que absorve luz, pero 
    // también la refleja más, pq la luz no puede penetrar y rebota, lo cual es clave para los colores metalicos. sin kapa
    // los materiales serían como vidrio. Sin eta perderían intensidad angular.
}
// --------------------------- HORA 8 ------------------------------------------------------------------------

//esta función es para generar las normales de microfacetas y sigue la distribución de Beckmann
//pra simular las superficies rugosas
//
Vector MuestreoMicrofacet(double aspereza, double u1, double u2) { // recibe aspereza, y dos números aleatorios
    double alpha = aspereza; // definimos la aspereza
    double tan2Theta = -alpha * alpha * log(1.0 - u1);  //  D(wh) = χ⁺(cos θ) * exp(-tan^2θ/α²) / (π α^2 cos^2θ)
                                                        // en la fórmula es arctan() pero si lo pasamos tenemos la tangente del ángulo polar
    double cos2Theta = 1.0 / (1.0 + tan2Theta); // 1 + tan^2θ = sec^2θ = 1/cos^2θ usando esta identidaad despejamos cos^2
    double cosTheta = sqrt(cos2Theta); // y sacamos la reíz cuadrada para tenerla normal

    double sin2Theta = 1.0 - cos2Theta; // sin^2θ + cos^2θ = 1  =>  sin^2θ = 1 - cos^2θ
    double sinTheta = sqrt(sin2Theta); // sacamos la raíz cuadrada para tenerla normal
    double phi = 2.0 * M_PI * u2; // y cálculamos el ángulo azimutal.

    return esfericasACartesianas(cosTheta, sinTheta, phi); // pasamos las coordenadas es esfpericas a cartecianas
    // y así obtenemos el vector wh
}

// en la 398 se calcula la reflactancia de un material conductor basado en el modelo de microfacetas
// aquí ya se comnina la D, G y F, pero sólo para encontrar la BRDF. Recibimos una esfera, el vector wi, wo y un vector n
Color BRDFMicrofacet(const Sphere &obj, const Vector &wi, const Vector &wo, const Vector &n) {
    Vector s, t, normal = n; // creamos el sistema de coordenasdas tomando el vector n 
    coordinateSystem(normal, s, t);  // tenemos que ponerlas en un sistema de coordenadas locales, ya que la 
                                    // distribución Beckmann asume que las microfacetas están orientadas a z como normal
                                    // el cálcullo par a las coordenadas locales es más simple. 
                                    // las fórmulas de microfacetas se escriben como si la superficie estuviera alineada con el eje Z
    
    Vector wiLocal = globalALocal(wi, normal, s, t); // hacemos wi local
    Vector woLocal = globalALocal(wo, normal, s, t); // y wo local
    
    if (wiLocal.z <= 1e-6 || woLocal.z <= 1e-6) return Color(); // asefuramos que wi y wo estén en eñ hemisferio correcto
                                                              // si es muy pequeño, entonces lo ponemos como sin color.
                                                              // en caso de alpha algo algunas wh podrían no ser validas
                                                              // aparte de que ambas se usarán como denominador y no puende ser 0
                                                              // si las dejamos como 0 podríamos tener artefactos negros o brillos incorrectos
    Vector whLocal = (wiLocal + woLocal).normalize(); // se saca wh local y se normaliza para sólo tener la dirección
                                                      
    
    double D = DistribucionBeckmann(whLocal, obj.aspereza); //sacamos D con wh local y la asperza
    double G = TerminoG(woLocal, wiLocal, whLocal, obj.aspereza); //también G
    Color F = FresnelConductor(abs(wiLocal.dot(whLocal)), obj.eta, obj.kappa); // y F se usa abs para evitar errores numéricos
    
    double denominador = 4.0 * abs(woLocal.z) * abs(wiLocal.z);  // se usa abs para evitar erroes
    if (denominador < 1e-6) return Color(); // para evitar dividir entre 0
    
    return Color(F.x * D * G / denominador, F.y * D * G / denominador, F.z * D * G / denominador);
} //Muestreo: "¿Qué direcciones son físicamente posibles?"
  //BRDF: "¿Cuánta luz viaja en esas direcciones?"
// --------------------------- HORA 9 ------------------------------------------------------------------------

// función general de Monte Carlo que acepta cualquier estrategia de muestreo
// monteCarlo acepta una estrategia de muestreo como parámetro, la cual es un apuntador a una función que
// acepta cualquier función que regrese un MuestreoResult, y tenga como parámetros un Point, Vector constantes y double
Color monteCarlo(int N, Vector &n, Point &x, const Sphere &obj,
                 MuestreoResult (*estrategiaMuestreo)(const Point&, const Vector&, double&)) {
    Color sum; // declaramos una variable para la suma
    Point puntoSeguro = x + n * 1e-4;  // se le agreg 1e-4 para evitar que haya intersecciones consigo misma por errors de cálculo
                                       // así que lo desplazamos un poco, generalmente fuera de la superficie.
                                       // al multiplicar después debemos ajustarlo con la esfera 
                                       // se le suma la normal para que el rayo quede afuera, ya que n tiene una dirección
                                       // hacia afuera de la esfera 
    for (int i = 0; i < N; ++i) {  // se usa ++i en lugar de i++ ya que i++ genera una copia innecesaria pq devuelve el valor original de 
                                   //i y luego lo incrementa,  y ++i incrementa su valaor y 
                                   // devuelve el valor utilizado. 
        double prob; // declaramos una variable de probabilidad
        MuestreoResult resultado = estrategiaMuestreo(puntoSeguro, n, prob); // se usa la estrategia de muestreo pasada y se guarda el resultado
                                                                             // el cual es la estructura del resultado
        Vector wi = resultado.wi; // optenemos el vector wi del resultado en otra variable para  que sea más legible

        Ray shadowRay(puntoSeguro, wi); // y lanzamos un rayo desde el punto seguro con dirección a wi
        double t_shadow; // declaramos variables para t_min
        int id_shadow; // y id_min

        bool hayInterseccion = intersect(shadowRay, t_shadow, id_shadow); // verificamos si sí hay al menos una intersección
                                                                          // en la escena con ere rayo. Si sí, es true y tenemos 
                                                                          // la t más cercana y el id de la esfera más cercana          

        if (!hayInterseccion) { // si ni hay intersección, pasamos a la siguiente muestra
            continue;
        }
        
        const Sphere &objIntersectado = spheres[id_shadow]; // si sí hubo, guardamos la esfera que se intersectó.
                                                            // se usa & para trabajar con la esferea original.
                                                            // es como un alias para esa esfera en específico
                                                            // no se le puede cambiar el valor.
        Color Le = objIntersectado.e; // guardamos la emición de ese objeto
        Color fr = obj.c * (1.0 / M_PI); // y cáculamos la BRDF la cual es lambertiana
        double cosTheta = n.dot(wi); // cálculamos coTheta para la fórmula

        if (cosTheta > 0) { // si cosTheta es menor, significa que no hay contribución por estar atrás
            Color contribucion = Le.mult(fr) * (cosTheta / prob); // si sí hay, el color de la contribucuón en ese pixel
                                                                  // es la multiplicación de Le*fr * cosTheta/prob
            sum = sum + contribucion; // entonces lo sumamos con lo que ya teníamos
        }
    }
    return sum * (1.0 / N); // y regresamos la suma entre la cantidad de muestras
}
// --------------------------- HORA 10 ------------------------------------------------------------------------
// al pasar el nombre de una fucnión como parámetro en automático se convierte en un puntero a una función
Color monteCarloHemisfer(int N, Vector &n, Point &x, const Sphere &obj) { // pide un número de muestras, un vector n, punto x y una esfera
    return monteCarlo(N, n, x, obj, MuestreoUniformeHemisferico); // se hace el montecarlo psándile el tipo de muestreo como parámetro
}

// igual pedimos la cantidad de muestras, un vector n, y una esfera y a monteCarlo le pasamos estos parámetros con el tipo de muestreo
Color monteCarloCoseno(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoCosenoHemisferico);
}

// las funciones wrapper son la envoltura para poder llamar a monteCarlo. 
MuestreoResult MuestreoAreaWrapper(const Point &x, const Vector &n, double &prob) {
    return MuestreoArea(x, prob); // el muestreo de área sólo ocupa dos argumentos ya que no ocupa buscar la fuente de luz
}

// est afucnión usa a monteCarlo para que llame a la envoltura, ya que se necesita el parámetro del vector n
// que el muestreo de área no necesita, entonces wrapper sí tiene n aunque no la use sólo para después llamar a
// la funci´n del muestreo
Color monteCarloArea(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoAreaWrapper);
}
// monteCarloAnguloSolido llama a monteCarlo, la cual llama a esta fución que es una envoltura para el muestreo en
// ángulo sólido, ya que en el puntero se piden argumentos especídifos y el ángulo sólido sólo ocupa dos de ellos
// entonces se pone el vector n en la envoltura para que monteCarlo acepte esta función como argumento
MuestreoResult MuestreoAnguloSolidoWrapper(const Point &x, const Vector &n, double &prob) {
    return MuestreoAnguloSolido(x, prob);
}

// llama a monteCarlo pasándole como argumento la envoltura de la función a usar
Color monteCarloAnguloSolido(int N, Vector &n, Point &x, const Sphere &obj) {
    return monteCarlo(N, n, x, obj, MuestreoAnguloSolidoWrapper);
}
// --------------------------- HORA 11 ------------------------------------------------------------------------

// el muestreo de fuente puntual en 477 se le pasa un vectior n, un punto y una esfera. 
// se hizo una función aparte de monteCarlo ya que para una fuente puntual hacemos el muestreo directo de la fuente de luz
// ya que sólo al haber un punto que emite luz, sólo hay una dirección de donde puede tener iluminación, así que 
// que se quita el bucle para la sumatoria.
Color monteCarloFuentePuntual(Vector &n, Point &x, const Sphere &obj) {
    Point puntoSeguro = x + n * 1e-4; // el punto es seguro ya que si lo dejamos tal cual, la fórmula lo toma como
                                      // un borde de la esfera perfecto, pero puede que esté un poquito más abajo y tenga una 
                                      // auto-intersección, entonces lo subimos 1e-4 en la dirección de la normal (arriba) para evitar 
                                      // que ocurra esto.
    // obtener dirección y punto
    double prob; // declaramos una variable para la probabilidad
    MuestreoResult resultado = MuestreoFuentePuntual(puntoSeguro, prob); // usamos el muestreo de la funte puntual
    Vector wi = resultado.wi; //y del resultado sacamos la wi calculada
    Point puntoEnFuente = resultado.puntoEnFuente; // declaramos un punto como el punto en fuente.
                                                   // estas las declaramos como variables aparte para no tener que 
                                                   // estar haciendo resultado.puntoEnFuente y sea más legible
    // calcular distancia a la fuente
    Vector direccionCompleta = wi; // sacamos otra variable para poder tener la dirección completa aunque está de más. se puede cambiar
    double distanciaAlCuadrado = direccionCompleta.dot(direccionCompleta); // sacamos la distancia al cuadrado ya que el producto punto nos da un escalar
    double distancia = sqrt(distanciaAlCuadrado); // y le sacamos raís para tener la distancia 
    wi = wi.normalize(); // y una vez hecho eso, normalizamos nuestro vector
    // verificar si hay obstáculos 
    Ray shadowRay(puntoSeguro, wi); // creamos un rato para ver si hay visibilidad
    double t_shadow; 
    int id_shadow;
    
    if (intersect(shadowRay, t_shadow, id_shadow)) { // si logra hacer una intersección en el camino hacia la luz
        // y esa intersección es menor a la distancia que hay hasta la luz
        if (t_shadow < distancia - 1e-4) {
            return Color(); // significa que no llegó y que hay algo tapandi la luz
        }
    }

    Color Le = fuenteLuminosa.e; // guardamos la emisión de la fuente luminosa
    Color fr = obj.c * (1.0 / M_PI); // BRDF lambertiana
    double cosTheta = n.dot(wi); // tenemos el cactor cosTheta
    
    if (cosTheta > 0) { // si este es mayor a 0
        // tomamos la contribución de la luz en ese punto, el cual es toda la que tendrá 
        // para fuente puntual, la intensidad decrece con 1/r^2
        Color contribucion = Le.mult(fr) * (cosTheta / distanciaAlCuadrado);
        return contribucion; // y regresamos la contribución 
    }
    
    return Color(); // no hay contribución si cosTheta <= 0
}
// --------------------------- HORA 12 ------------------------------------------------------------------------


// es path tracing explícito, ya que se toma la iluminación directa de manera explícita y apartir de allí genera un nuevi rayo según 
// el material. Primero verificamos si hay intersecciones, vemos el tipo de muestreo para iluminación directa, usamos el tipo
// de muestreo según el material, y usamos recursión para generar nuevos rayos 
Color pathTracing(const Ray &ray, int profundidad, int maxProfundidad) { // recibimos un rayo, la profunfidad donde nos encontramos y la profundidad máxima
    if (profundidad >= maxProfundidad) return Color(); // si se llega al punto en el que la profundidad 
                                                       // actual ya llegó a la máxima permitida, se regresa un color negro.
    double t; // definimos t para las intersecciones
    int id; // e id para la esfera más cercana
    if (!intersect(ray, t, id)) return Color(); // si no hay intersección, se regresa un color negro (para nuestra escena no puede pasar ya que está rodeada en una cajita)
    
    const Sphere &obj = spheres[id]; // guardamos la esfera más cercana 
    Point x = ray.o + ray.d * t; // y cálculamos el punto de la intersección, el cuál es el origen del rayo, con la 
                                 // dirección del rayo por la distancia hacia la esfera
    Vector n = (x - obj.p).normalize(); // cálculamos la normal que va desde el centro de la esfera al punto de intersección y la normalizamos
    
    // Si es fuente luminosa, retornar emisión
    if (obj.e.x > 0 || obj.e.y > 0 || obj.e.z > 0) { // se verifican todos los valores de la emisión ya que 
        return obj.e;                                // podría pasar que se  emita luz de otros colores
    }
 
    Point puntoSeguro = x + n * 1e-4; // el punto seguro lo sacamos para evitar autointersección, ya que matemáticamente es una esfera
                                      // con un borde perfecto, pero si el punto está un poquito más abajo se puede tomar
                                      // como una falta intersección, aí que lo subimos tantito
    Color L(0,0,0); // definimos un color L el cual será la emisión total de ese punto
    
    // Iluminación directa. Se puede cambiar el tipo de muestreo
    L = L + monteCarloAnguloSolido(N_,n, x, obj);  // el el código tenemos varios tipos de muestreo, si queremos usar path tracing para
                                             // cualquier tipo de muestreo, se debe de cambiar aquí. Tendríamos muestreo de área
                                             // de coseno, fuente puntual, etc
    Vector wo = Vector(0,0,0) - ray.d; // declaramos un vector de entrada, el cual va desde 0 con la dirección del rayo, pero invversa
                                       // wo es la dirección desde el punto de intersección hasta el observador. Es la dirección por la cual la luz "sale" del punto hacia donde estamos mirando.
                                       // al restarle a un vector 0 la dirección del rayo, hacemos que la dirección esté al lado contrario
                                       // ahora sí saliendo desde la intersección hasta la camara   
    if (obj.material == DIFFUSO) { // si el material de la esferea es difuso
        double prob; // declaramos una probabilidad
        MuestreoResult resultado = MuestreoCosenoHemisferico(puntoSeguro, n, prob); // usamos el muestreo de conseno hemisferico ya que es el menor para materiales difusos
                                                                                    //n  es el vector que sale de la primera intersección
        Vector wi = resultado.wi; // y optenemos del resultado la dirección wi
// --------------------------- HORA 12 ------------------------------------------------------------------------
    
        if (prob > 1e-6 && n.dot(wi) > 1e-6) {  // verificamos si la probabilidad es más grande que 0, ya que en el muestreo de área
                                                // o en ángulo sólido podría haber estos casos
                                                // tambien se garantiza que wi no sea perpendicular a n, esto puede pasar en el muestreo 
                                                // de coseno hemisferico
            Ray nuevoRayo(puntoSeguro, wi);  // creamos un nuevo rayo con la nueva dirección que obtuvimos con el muestreo
            Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad); // hacemos las llamadas recursivas
                                                                                // para los rebotes, sumando uno al reppte actual
            Color fr = obj.c * (1.0 / M_PI); // sacamos la BRDF para el material lambertiano
            double cosTheta = n.dot(wi); // sacamos el coseno de theta para esa dirección y la normal
            
            double factor = (cosTheta / prob); // es el factor para corregir el sesgo introducido por el muestreo
            L = L + Color(fr.x * Li.x * factor, fr.y * Li.y * factor, fr.z * Li.z * factor); //L es el nuevo valor de la radiancia, que contiene
                                                                            // la iluminación directa más el color dado por la iluminación global
                                                                            // el cual es el valor de la BRDF por la ilumincación total
                                                                            // de los caminos por el factor del sesgo 
        }
    } 
    else if (obj.material == CONDUCTOR) { // si el material es un conductor
        double u1 = (double)rand() / RAND_MAX; // definimos números aleatorios 
        double u2 = (double)rand() / RAND_MAX; 
        
        Vector s, t;
        coordinateSystem(n, s, t); // creamos un sistema de coordenadas con respecto a n
        
        Vector whLocal = MuestreoMicrofacet(obj.aspereza, u1, u2); // se hace el muestreo de microfacet para tener 
                                                                   // wh de manera local
        Vector wh = localAGlobal(whLocal, n, s, t); // ahora este vector lo pasamos a global
        wh = wh.normalize(); // se normaliza sólo para tener la dirección
        
        Vector wi = wh * (2.0 * wo.dot(wh)) - wo; // es la relexión especular sobre la micronormal de una faceta
                                                  // que se usa para modelar materiales conductores 
                                                  // con esto nos aseguramos que el vector wi tenga el mismo ángulo pero contrario que wo
                                                  // wo.dot(wh) es el coseno del ángulo entre wo y wh
                                                  //  wo.dot(wh) * wh es la proyección de wo en wh. La fórmula se despeja encontrando
                                                  // un vector perpendicular a wo.
        wi = wi.normalize(); // lo normalizamos ya que sólo queremos la dirección
        
        if (n.dot(wi) > 1e-6) { // si no nos aseguramos que wi sí esté del otro lado de n, entonces puede dar ruido
            double D = DistribucionBeckmann(whLocal, obj.aspereza); // cálculamos cuántas microfacetas hay alineadas a n
            double G = TerminoG(globalALocal(wo, n, s, t), globalALocal(wi, n, s, t), whLocal, obj.aspereza); // vemos cuáles de esas sí contribuyen 
            Color F = FresnelConductor(fabs(wo.dot(wh)), obj.eta, obj.kappa); // y vemos la reflectancia. 
                                                                              // se usa fabs par garantizar que sea positivo
            
            // Cálculo de la probabilidad de muestreo de la dirección ωₒ dado ωₕ (microfaceta)
            // - D: Valor de la NDF (Normal Distribution Function) en ωₕ (densidad de microfacetas)
            // - whLocal.z: Componente z de ωₕ en coordenadas locales (equivale a n·ωₕ)
            // - wo.dot(wh): Coseno del ángulo entre ωₒ y ωₕ (para el Jacobiano 4|ωₒ·ωₕ|)
            // La fórmula es: p(ωₒ) = D(ωₕ) * (n·ωₕ) / (4|ωₒ·ωₕ|) (Tema 5, p. 48)
            double probMicrofacet = D * whLocal.z / (4.0 * fabs(wo.dot(wh)));

            // Denominador ajustado para el estimador Monte Carlo:
            // - probMicrofacet: p(ωₒ) calculada arriba
            // - n.dot(wh): Producto punto entre la normal macroscópica (n) y ωₕ (microfaceta)
            // Se multiplica p(ωₒ) * (n·ωₕ) para cancelar términos en la BRDF posteriormente
            // Nota: fabs() garantiza no división por negativos (conservación de energía)
            double denom = probMicrofacet * fabs(n.dot(wh));

            // Evitar división por cero o valores numéricamente inestables (< 1e-6)
            if (denom > 1e-6) { 
                // Generar nuevo rayo para path tracing recursivo:
                // - puntoSeguro: Punto de intersección desplazado en la dirección de n (para evitar self-intersection)
                // - wi: Dirección muestreada (ωᵢ = 2(ωₕ·ωₒ)ωₕ - ωₒ, reflexión especular sobre ωₕ)
                Ray nuevoRayo(puntoSeguro, wi);

                // Llamada recursiva para calcular la radiancia entrante (Lᵢ):
                // - nuevoRayo: Rayo hacia la dirección muestreada wi
                // - profundidad + 1: Incrementar contador de rebotes
                // - maxProfundidad: Límite de rebotes para evitar recursión infinita
                Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad);

                // Factor de corrección Monte Carlo:
                // - G: Término geométrico (masking-shadowing) de la BRDF
                // - wo.dot(wh): |ωₒ·ωₕ| (aparece en el numerador de la BRDF)
                // - denom: p(ωₒ)*(n·ωₕ) (denominador ajustado)
                // Fórmula final: (G * |ωₒ·ωₕ|) / (p(ωₒ)*(n·ωₕ)) = (G * 4|ωₒ·ωₕ|²) / (D*(n·ωₕ)²) (tras sustituir p(ωₒ))
                double factor = (G * fabs(wo.dot(wh)) / denom);

                // Acumular contribución de la muestra actual a la radiancia total (Lₒ):
                // - F: Término de Fresnel (dependiente del material)
                // - Li: Radiancia entrante calculada recursivamente
                // - factor: Corrección de importancia muestral
                // Se descompone por canales RGB (x, y, z) para preservar el color
                L = L + Color(F.x * Li.x * factor, F.y * Li.y * factor, F.z * Li.z * factor);
            }
        }
    }
    
    return Color(fmin(fmax(L.x, 0.0), 1.0), fmin(fmax(L.y, 0.0), 1.0), fmin(fmax(L.z, 0.0), 1.0));
}

// es una función general que nos dice si hubi una intersección es una esfera con un rayo
// y guarda el punto de intersección, y la normal en el punto de intersección
bool shade(const Ray &r, Point &x, Vector &n, const Sphere *&obj) {  // recibe un rayo, un punto, un vector n y una esfera
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

// es una función que usa shade para el muestreo de hemisferio
Color shadeHemisfer(const Ray &r) { // recibe un rayo
    Point x; // de claramos un punto x
    Vector n; // un vector n
    const Sphere *obj; // y una esfera
    
    if (!shade(r, x, n, obj)) // vemos si hay 
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