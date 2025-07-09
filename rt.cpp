// para cambiar de muestreo es con las funciones shade
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  
#include <omp.h>
#include <time.h>


#define N_ 32
#define N_ESFERAS 8
#define N_ARREGLO 7

class Vector 
{
public:        
	double x, y, z; // coordenadas x,y,z 
  
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
typedef Vector Point;
typedef Vector Color;

enum MaterialType {
    DIFFUSO,
    CONDUCTOR
};

class Ray 
{ 
public:
	Point o;
	Vector d; // origen y direcccion del rayo
	Ray(Point o_, Vector d_) : o(o_), d(d_) {} // constructor
};
 
class Sphere 
{
public:
	double r;	// radio de la esfera
	Point p;	// posicion
	Color c;	// color  
	Color e;
    MaterialType material;
    double aspereza; // para microfacet
    Color eta, kappa; // para conductor

	Sphere(double r_, Point p_, Color c_, Color e_ = Color(), MaterialType mat = DIFFUSO, double aspe = 0.1, Color eta_ = Color(1,1,1), Color kappa_ = Color(0,0,0)) 
        : r(r_), p(p_), c(c_), e(e_), material(mat), aspereza(aspe), eta(eta_), kappa(kappa_) {}
  
	double intersect(const Ray &ray) const {
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

class FuenteLuminosa : public Sphere {
public: 
	FuenteLuminosa(double r_, Point p_, Color c_, Color e_ = Color(10, 10,10))
		: Sphere(r_, p_, c_, e_) {}	
};

class FuentePuntual : public Sphere {
public:
    FuentePuntual(Point p_, Color e_ = Color(4000, 4000, 4000), double r_ = 0) 
        : Sphere(r_, p_, Color(), e_) {} 
};

struct MuestreoResult {
	Vector wi; 
	double prob;
	Point puntoEnFuente;

	MuestreoResult(Vector w, double p, Point pt = Point()) : wi(w), prob(p), puntoEnFuente(pt) {}
};

Sphere spheres[] = {
    //Escena: radio, posicion, color, emisión
     	Sphere(1e5,  Point(-1e5 - 49, 0, 0),     Color(.75, .25, .25)), // pared izq
        Sphere(1e5,  Point(1e5 + 49, 0, 0),     Color(.25, .25, .75)), // pared der
        Sphere(1e5,  Point(0, 0, -1e5 - 81.6),  Color(.25, .75, .25)), // pared detras
        Sphere(1e5,  Point(0, -1e5 - 40.8, 0),  Color(.25, .75, .75)), // suelo
        Sphere(1e5,  Point(0, 1e5 + 40.8, 0),  Color(.75, .75, .25)), // techo
        // esferas normales para el resto de muestreos
        /*Sphere(16.5, Point(-23, -24.3, -34.6),  Color(.2, .3, .4)), 
        Sphere(16.5, Point(23, -24.3, -3.6),     Color(.4, .3, .2)),*/

        // esferas de aluminio y oro para microfacet
        Sphere(16.5, Point(-23, -24.3, -34.6),  Color(.2, .3, .4), Color(), CONDUCTOR, 0.3, 
           Color(1.44, 0.96, 0.61), Color(7.47, 6.52, 5.29)), // Aluminio
        Sphere(16.5, Point(23, -24.3, -3.6),     Color(.4, .3, .2), Color(), CONDUCTOR, 0.3,
           Color(0.143, 0.374, 1.442), Color(3.982, 2.386, 1.603)), // Oro */

        // fuentes luminosas   
        FuenteLuminosa(10.5, Point(0, 24.3, 0), Color(1, 1, 1)) // fuente normal
    	//FuentePuntual(Point(0, 24.3, 0)) // fuente puntual 

       // FuenteLuminosa(10.5, Point(-23, 24.3, 0), Color(1, 1, 1), Color(12, 5, 5)),
       // FuenteLuminosa(5, Point(23, 24.3, -50), Color(1, 1, 1), Color(5, 5, 12)) 
};

Sphere& fuenteLuminosa = spheres[N_ARREGLO];

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

void coordinateSystem(Vector &n, Vector &s, Vector &t){
	if (fabs(n.x) > fabs(n.y))
	{
		float invLen = 1.0f / sqrt(n.x * n.x + n.z * n.z);
		t = Vector(n.z * invLen, 0.0f, -n.x * invLen);
	} else {
		float invLen = 1.0f / sqrt(n.y * n.y + n.z * n.z);
		t = Vector(0.0f, n.z * invLen, -n.y * invLen);
	}
	s = t%n;
}

Vector esfericasACartesianas(double cosTheta, double sinTheta, double phi){
	return Vector(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}

Vector localAGlobal(const Vector &local, Vector &n, Vector &s, Vector &t){
	return Vector(
		s.x * local.x + t.x * local.y + n.x * local.z,
		s.y * local.x + t.y * local.y + n.y * local.z,
		s.z * local.x + t.z * local.y + n.z * local.z 
	);
}

Vector globalALocal(const Vector &global, Vector &n, Vector &s, Vector &t){
	return Vector(
		global.dot(s),
		global.dot(t),
		global.dot(n)
	);
}

MuestreoResult MuestreoUniformeHemisferico(const Point &x, const Vector &n, double &prob) {
	double u1 = (double)rand() / RAND_MAX;
	double u2 = (double)rand() / RAND_MAX;

	double theta = acos(u1); // formulas del muestreo
	double phi = 2.0 * M_PI * u2;

	double cosTheta = cos(theta); // para las coordenadas locales 
	double sinTheta = sin(theta);

	Vector s, t, normal = n; // delclaramos los vectores para hacer el sistema de coordenadas local
	coordinateSystem(normal, s, t);

	Vector localDir = esfericasACartesianas(cosTheta, sinTheta, phi); //las coordenadas que tenemos son esféricas
	Vector wi = localAGlobal(localDir, normal, s, t); // las paramos a global
	wi = wi.normalize(); //y lo normalizamos por que sólo nos interesa la dirección
	 
	prob = 1.0 / (2.0 * M_PI); // probabilidad uniforme sobre el hemisferio

	return MuestreoResult(wi, prob, Point());
}


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

double xMas(double x) {
    return x > 0 ? 1.0 : 0.0;
}

double DistribucionBeckmann(const Vector &wh, double aspereza) {
    double cosTheta = wh.z; // ángulo con la normal
    double alpha = aspereza;
    double cos2Theta = cosTheta * cosTheta;
    double tan2Theta = (1.0 - cos2Theta) / cos2Theta;
    
    double x = xMas(cosTheta);
    return x * exp(-tan2Theta / (alpha * alpha)) / 
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
    Color Rpe = Color((t1.x - t2.x) / (t1.x + t2.x),
                    (t1.y - t2.y) / (t1.y + t2.y),
                    (t1.z - t2.z) / (t1.z + t2.z));
    
    Color t3 = Color(cos2Theta * a2Masb2.x + sin2Theta * sin2Theta,
                    cos2Theta * a2Masb2.y + sin2Theta * sin2Theta,
                    cos2Theta * a2Masb2.z + sin2Theta * sin2Theta);
    
    Color t4 = t2 * sin2Theta;
    Color Rpa = Color(Rpe.x * (t3.x - t4.x) / (t3.x + t4.x),
                    Rpe.y * (t3.y - t4.y) / (t3.y + t4.y),
                    Rpe.z * (t3.z - t4.z) / (t3.z + t4.z));
    
    return Color(0.5 * (Rpe.x + Rpa.x), 0.5 * (Rpe.y + Rpa.y), 0.5 * (Rpe.z + Rpa.z));
}

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
    L = L + monteCarloAnguloSolido(N_,n, x, obj);
    
    // Muestreo según material
    Vector wo = Vector(0,0,0) - ray.d;
    
    if (obj.material == DIFFUSO) {
        double prob;
        MuestreoResult resultado = MuestreoCosenoHemisferico(puntoSeguro, n, prob);
        Vector wi = resultado.wi;
        
        if (prob > 1e-6 && n.dot(wi) > 1e-6) {
            Ray nuevoRayo(puntoSeguro, wi);
            Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad);
            Color fr = obj.c * (1.0 / M_PI);
            double cosTheta = n.dot(wi);
            
            double factor = (cosTheta / prob);
            L = L + Color(fr.x * Li.x, fr.y * Li.y, fr.z * Li.z) * factor;
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
            
            if (denom > 1e-6) {
                Ray nuevoRayo(puntoSeguro, wi);
                Color Li = pathTracing(nuevoRayo, profundidad + 1, maxProfundidad);
                double factor = (G * fabs(wo.dot(wh)) / denom);
                
                L = L + Color(F.x * Li.x, F.y * Li.y, F.z * Li.z) * factor;
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

Color shadeMicrofacet(const Ray &r) {
    Color L(0,0,0);
    int spp = 64; // samples por pixel
    for (int i = 0; i < spp; ++i) {
        Color sample = pathTracing(r, 0, 10);
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