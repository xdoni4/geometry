#include <iostream>
#include <math.h>
#include <vector>
#include <initializer_list>
#include <algorithm>

struct Constants {
    static constexpr double eps = 1e-3;
    static constexpr double pi = 3.1415926;
};

struct Point {
    double x, y;
    Point(): x(0), y(0){}
    Point(double x, double y): x(x), y(y){}
    Point (const Point &p2): x(p2.x), y(p2.y){}
    Point& operator=(const Point &p2) {
        x = p2.x;
        y = p2.y;
        return (*this);
    }

    Point& operator+=(const Point &p2) {
        x += p2.x;
        y += p2.y;
        return (*this);
    }

    Point& operator-=(const Point &p2) {
        x -= p2.x;
        y -= p2.y;
        return (*this);
    }

    Point& operator*=(double k) {
        x *= k;
        y *= k;
        return (*this);
    }

    Point& operator/=(double k) {
        x /= k;
        y /= k;
        return (*this);
    }
};

Point operator+(const Point &p1, const Point &p2) {
    Point p(p1);
    p += p2;
    return p;
}

Point operator-(const Point &p1, const Point &p2) {
    Point p(p1);
    p -= p2;
    return p;
}

Point operator*(const Point &p1, double k) {
    Point p(p1);
    p *= k;
    return p;
}

Point operator*(double k, const Point &p1) {
    Point p(p1);
    p *= k;
    return p;
}

Point operator/(const Point& p1, double k) {
    Point p(p1);
    p /= k;
    return p;
}

double dist(const Point &p1, const Point &p2) {
    Point p = p1 - p2;
    return sqrt(p.x * p.x + p.y * p.y);
}

bool operator==(const Point &a, const Point &b) {
    return (fabs(a.x - b.x) < Constants::eps) && (fabs(a.y - b.y) < Constants::eps);
}

bool operator!=(const Point &a, const Point &b) {
    return !(a == b);
}

std::ostream& operator<<(std::ostream &out, const Point &p) {
    out << '(' << p.x << ", " << p.y << ')'; 
    return out;
}

class Line {
public:
    double k, b, x;
    Line(const Point &p1, const Point &p2) {
        k = p2.x > p1.x ? (p2.y - p1.y) / (p2.x - p1.x) : (p1.y - p2.y) / (p1.x - p2.x);
        b = (p2.y - k * p2.x);
        x = p1.x;
    }

    Line(double k, double b): k(k), b(b) {}

    Line(const Point &p, double koef) {
        k = koef;
        b = (p.y - k * p.x);
        x = p.x;
    }
};

Point lineIntersection(const Line &l1, const Line &l2) {
    if (fabs(l1.k) > 1/Constants::eps && fabs(l2.k) < 1/Constants::eps) return Point(l1.x, l2.k*l1.x+l2.b);
    else if (fabs(l2.k) > 1/Constants::eps && fabs(l1.k) < 1/Constants::eps) return Point(l2.x, l1.k*l2.x+l1.b);
    double x = (l2.b - l1.b) / (l1.k - l2.k);
    double y = l1.k * x + l1.b;
    return Point(x, y);
}

bool operator==(const Line &l1, const Line &l2) {
    if (fabs(l1.k) > 1 / Constants::eps && fabs(l2.k) > 1 / Constants::eps) {
        if (l1.x == l2.x) return true;
        return false;
    }
    return (fabs(l1.k - l2.k) < Constants::eps) && (fabs(l1.b - l2.b) < Constants::eps);
}

bool operator!=(const Line &l1, const Line &l2) {
    return !(l1 == l2); 
}

std::ostream& operator<<(std::ostream &out, const Line &l) {
    out << l.k << "x + " << l.b;
    return out;
}

class Shape {
public:
    double crossProduct(const Point &p1, const Point &p2) const {
        return p1.x * p2.y - p1.y * p2.x;
    }

    double scalarProduct(const Point &p1, const Point &p2) const {
        return p1.x * p2.x + p1.y * p2.y;
    }

    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape &another) const = 0;
    virtual bool operator!=(const Shape &another) const = 0;
    virtual bool isCongruentTo(const Shape &another) const = 0;
    virtual bool isSimilarTo(const Shape &another) const = 0;
    virtual bool containsPoint(const Point &point) const = 0;
    virtual bool rotate(const Point &center, double angle) = 0;
    virtual bool reflect(const Point &center) = 0;
    virtual bool reflect(const Line &axis) = 0;
    virtual bool scale(const Point &center, double coefficient) = 0;
    virtual ~Shape() {}
};

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;

public:
    Polygon(const std::vector<Point> &v) { 
        vertices = v;
        double miny = 1e10;
        int b = -1;
        for (size_t i = 0; i < vertices.size(); i++) {
            if (vertices[i].y < miny) {
                miny = vertices[i].y;
                b = i;
            }
        }
        if (vertices.size() >= 3) {
            if (crossProduct(vertices[(b+1)%vertices.size()] - vertices[b], vertices[(vertices.size()+b-1)%vertices.size()] - vertices[b]) < -Constants::eps) std::reverse(vertices.begin(), vertices.end());
        }
    }

    Polygon(const std::initializer_list<Point> &lst) {
        vertices = lst;
        double miny = 1e10;
        int b = -1;
        for (size_t i = 0; i < vertices.size(); i++) {
            if (vertices[i].y < miny) {
                miny = vertices[i].y;
                b = i;
            }
        }
        if (vertices.size() >= 3) {
            if (crossProduct(vertices[(b+1)%vertices.size()] - vertices[b], vertices[(vertices.size()+b-1)%vertices.size()] - vertices[b]) < -Constants::eps) std::reverse(vertices.begin(), vertices.end());
        }
    } 

    template<class ... Points>
    Polygon(const Points& ... points) : vertices{points ...} {
        double miny = 1e10;
        int b = -1;
        for (size_t i = 0; i < vertices.size(); i++) {
            if (vertices[i].y < miny) {
                miny = vertices[i].y;
                b = i;
            }
        }
        if (vertices.size() >= 3) {
            if (crossProduct(vertices[(b+1)%vertices.size()] - vertices[b], vertices[(vertices.size()+b-1)%vertices.size()] - vertices[b]) < -Constants::eps) std::reverse(vertices.begin(), vertices.end());
        }
    }

    std::vector<Point> getVertices() const {
        return vertices;
    }

    bool isConvex() const {
        bool sign = crossProduct(vertices[vertices.size()-1] - vertices[vertices.size()-2], vertices[0] - vertices[vertices.size()-1]) > 0;
        if (sign != (crossProduct(vertices[0] - vertices[vertices.size()-1], vertices[1] - vertices[0]) > 0)) return false;
        for (size_t i = 0; i < vertices.size() - 2; i++) {
            bool cursign = crossProduct(vertices[i+1] - vertices[i], vertices[i+2] - vertices[i+1]) > 0;
            if (cursign != sign) return false;
        }
        return true;
    }

    double perimeter() const override {
        double res = 0;
        for (size_t i = 0; i < vertices.size(); i++) {
            res += dist(vertices[i], vertices[(i+1)%vertices.size()]);
        }
        return res;
    }
    
    double angle(size_t i) const {
        return acos(scalarProduct(vertices[(i+1)%vertices.size()] - vertices[i], vertices[(vertices.size()+i-1)%vertices.size()] - vertices[i]) / dist(vertices[(i+1)%vertices.size()], vertices[i]) / dist(vertices[(vertices.size()+i-1)%vertices.size()], vertices[i]));
    }

    double area() const override {
        double res = 0;
        for (size_t i = 0; i < vertices.size(); i++) {
            res += crossProduct(vertices[i], vertices[(i+1)%vertices.size()]);
        }
        return fabs(res) / 2;
    }

    bool operator==(const Shape &another) const override {
        const Shape *ptr = &another;
        const Polygon *casted = dynamic_cast<const Polygon*>(ptr);
        if (casted == nullptr) return false;
        if (casted->vertices.size() != vertices.size()) return false; 

        size_t j = 0;
        size_t p = 0;
        for (size_t i = 0; i < vertices.size(); ) {
            while (vertices[i%vertices.size()] == casted->vertices[j%vertices.size()]) {
                i++; j++; p++;
                if (p == casted->vertices.size()) return true;
            }
            p = 0;
            j = 0;
            i++;
        }
        return false;
    }

    bool operator!=(const Shape &another) const override {
        return !((*this) == another);
    }

    bool isCongruentTo(const Shape &another) const override {
        const Shape *ptr = &another;
        const Polygon *casted = dynamic_cast<const Polygon*>(ptr);
        if (casted == nullptr) return false;
        if (vertices.size() != casted->vertices.size()) return false;
        
        size_t j = 0;
        size_t p = 0;
        int begin = -1;
        for (size_t i = 0; i < vertices.size(); ) {
            while (fabs(dist(vertices[i%vertices.size()], vertices[(i+1)%vertices.size()]) - dist(casted->vertices[j], casted->vertices[(j+1)%vertices.size()])) < Constants::eps) {
                i++; j++; p++;
                if (p == vertices.size()) {
                    begin = i%vertices.size();
                    break;
                }
            }
            p = 0;
            j = 0;
            i++;
        }
        if (begin == -1) return false;
        for (size_t i = 0; i < vertices.size(); i++) {
            if (fabs(angle((begin+i)%vertices.size()) - casted->angle(i)) >= Constants::eps) return false;
        } 

        return true;
    }

    bool isSimilarTo(const Shape &another) const override {
        const Shape *ptr = &another;
        const Polygon *casted = dynamic_cast<const Polygon*>(ptr);
        if (casted == nullptr) return false;
        if (vertices.size() != casted->vertices.size()) return false;
 
        size_t j = 0;
        size_t p = 0;
        int begin = -1;
        for (size_t i = 0; i < vertices.size(); ) {
            while (fabs(angle(i%vertices.size()) - casted->angle(j%vertices.size())) < Constants::eps) {
                i++; j++; p++;
                if (p == vertices.size()) {
                    begin = i%vertices.size();
                    break;
                }
            }
            p = 0;
            j = 0;
            i++;
        }
        if (begin == -1) return false;
        double _ratio = dist(casted->vertices[0], casted->vertices[1]) / dist(vertices[(begin+1)%vertices.size()], vertices[begin%vertices.size()]);  
        for (size_t i = 0; i < vertices.size(); i++) {
            if (fabs(_ratio*dist(vertices[(begin+i+1)%vertices.size()], vertices[(begin+i)%vertices.size()]) - dist(casted->vertices[i], casted->vertices[(i+1)%vertices.size()])) >= Constants::eps) return false;
        }
        return true;
    }

    bool containsPoint(const Point &point) const override { 
        Line l(point, 3.73);
        size_t intersections = 0, reps = 0;
        for (size_t i = 0; i < vertices.size(); i++) {
            Line curline(vertices[i], vertices[(i+1)%vertices.size()]);
            Point inter = lineIntersection(l, curline);
            double minx = std::min(vertices[i].x, vertices[(i+1)%vertices.size()].x), maxx = std::max(vertices[i].x, vertices[(i+1)%vertices.size()].x);
            double miny = std::min(vertices[i].y, vertices[(i+1)%vertices.size()].y), maxy = std::max(vertices[i].y, vertices[(i+1)%vertices.size()].y);
            bool a = (inter.x > minx-Constants::eps);
            bool b = (inter.x < maxx+Constants::eps);
            bool c = (inter.y > miny-Constants::eps);
            bool d = (inter.y < maxy+Constants::eps);
            bool e = (inter == vertices[i]);
            bool f = (inter == vertices[(i+1)%vertices.size()]);
            bool g = (inter.x > point.x);
            if (((a && b && c && d) || (e || f)) && g) intersections++;
            if ((inter == vertices[i]) || (inter == vertices[(i+1)%vertices.size()])) reps++;
        }
        intersections -= reps / 2; 
        return intersections & 1; 
    }
    
    bool rotate(const Point &center, double angle) override {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i] = Point(center.x + (vertices[i].x-center.x)*cos(angle*Constants::pi/180) - (vertices[i].y-center.y)*sin(angle*Constants::pi/180), center.y + (vertices[i].x-center.x)*sin(angle*Constants::pi/180) + (vertices[i].y-center.y)*cos(angle*Constants::pi/180));
        }
        return 0;
    }

    bool reflect(const Point &center) override {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i] = center - (vertices[i] - center);
        }
        return 0;
    }

    bool reflect(const Line &axis) override {
        for (size_t i = 0; i < vertices.size(); i++) {
            Point ort(axis.k, -1);
            if (vertices[i].y > axis.k*vertices[i].x + axis.b - Constants::eps)  ort *= -1;
            vertices[i] -= 2 * fabs(crossProduct(vertices[i] - Point(0, axis.b), Point(-1, -axis.k)))/(axis.k*axis.k+1) * ort;
        }
        return 0;
    }

    bool scale(const Point &center, double coefficient) override {
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] + (vertices[i] - center) * (coefficient - 1);
        }
        return 0;
    }
};

class Ellipse : public Shape {
protected:
    double e, a, b;
    Point f1, f2;

public:
    Ellipse(const Point &p1, const Point &p2, double r) {
        f1 = p1; 
        f2 = p2;
        e = dist(f1, f2) / r; 
        a = r / 2;
        double c = e*a;
        b = sqrt(a*a - c*c);
    }

    virtual ~Ellipse() {}

    std::pair<Point, Point> focuses() const {
        return std::make_pair(f1, f2);
    }

    std::pair<Line, Line> directrices() const {
        Point center = (f1 + f2) / 2; 
        Line l(f1, f2);
        double k_norm = -1/l.k;
        Point p1 = center + (f1-center)/(e*e);
        Point p2 = center + (f2-center)/(e*e);
        return std::make_pair(Line(p1, k_norm), Line(p2, k_norm));
    }

    double eccentricity() const {
        return e;
    }

    double perimeter() const override {
        return Constants::pi*(a + b)*(1 + (3 * (a-b)/(a+b) * (a-b)/(a+b) ) / (10 + sqrt(4 - (3 * (a-b)/(a+b) * (a-b)/(a+b) )))); 
    }

    double area() const override {
        return Constants::pi * a * b;
    }

    bool operator==(const Shape &another) const override {
        const Shape *ptr = &another;
        const Ellipse *casted = dynamic_cast<const Ellipse*>(ptr);
        if (casted == nullptr) return false;
        return f1 == casted->f1 && f2 == casted->f2 && fabs(a - casted->a) < Constants::eps && fabs(b - casted->b) < Constants::eps;
    }

    bool operator!=(const Shape &another) const override {
        return !((*this) == another);
    }

    bool isCongruentTo(const Shape &another) const override {
        const Shape *ptr = &another;
        const Ellipse *casted = dynamic_cast<const Ellipse*>(ptr);
        if (casted == nullptr) return false;
        return fabs(a - casted->a) < Constants::eps && (b - casted->b);
    }

    bool isSimilarTo(const Shape &another) const override {
        const Shape *ptr = &another;
        const Ellipse *casted = dynamic_cast<const Ellipse*>(ptr);
        if (casted == nullptr) return false;
        double _ratio = casted->a / a;
        return fabs(_ratio * b - casted->b) < Constants::eps ? true : false;
    }

    bool containsPoint(const Point &point) const override { 
        return (dist(f1, point) + dist(f2, point) < 2*a+Constants::eps);
    }

    bool rotate(const Point &center, double angle) override {
        f1 = Point(center.x + (f1.x-center.x)*cos(angle*Constants::pi/180) - (f1.y-center.y)*sin(angle*Constants::pi/180), center.y + (f1.x-center.x)*sin(angle*Constants::pi/180) + (f1.y-center.y)*cos(angle*Constants::pi/180));
        f2 = Point(center.x + (f2.x-center.x)*cos(angle*Constants::pi/180) - (f2.y-center.y)*sin(angle*Constants::pi/180), center.y + (f2.x-center.x)*sin(angle*Constants::pi/180) + (f2.y-center.y)*cos(angle*Constants::pi/180));
        return 0;
    }

    bool reflect(const Point &center) override {
        f1 = center - (f1 - center);
        f2 = center - (f2 - center);
        return 0;
    }

    bool reflect(const Line &axis) override {
        Point ort1(axis.k, -1);
        if (f1.y > axis.k*f1.x + axis.b - Constants::eps)  ort1 *= -1;
        f1 -= 2 * fabs(crossProduct(f1 - Point(0, axis.b), Point(-1, -axis.k)))/(axis.k*axis.k+1) * ort1;
        Point ort2(axis.k, -1);
        if (f2.y > axis.k*f2.x + axis.b - Constants::eps)  ort2 *= -1;
        f2 -= 2 * fabs(crossProduct(f2 - Point(0, axis.b), Point(-1, -axis.k)))/(axis.k*axis.k+1) * ort2;
        return 0;
    }

    bool scale(const Point &center, double coefficient) override {
        f1 = f1 + (f1 - center) * (coefficient - 1);
        f2 = f2 + (f2 - center) * (coefficient - 1);
        a *= coefficient;
        b *= coefficient;
        return 0;
    }
};

class Circle : public Ellipse {
protected:
    double r;
    Point c;

public:
    Circle(const Point &p1, double r) : Ellipse(p1, p1, 2*r), r(r), c(p1) {}
    virtual ~Circle() {}
    
    double radius() {
        return r;
    }

    Point center() {
        return c;
    }
};

class Rectangle : public Polygon {
public:
    Rectangle(const Point &p1, const Point &p2, double _ratio) : Polygon() {
        Line diag(p1, p2);
        double k = fabs(diag.k) > 1/Constants::eps ? -1/_ratio : (_ratio + diag.k) / (1 - _ratio * diag.k);
        Line l1(p1, k);
        Line l2(p2, -1/k);
        Point p3 = lineIntersection(l1, l2);
        Point p4 = p1 + (p2 - p3);
        vertices = {p1, p3, p2, p4};
    }

    virtual ~Rectangle() {}

    Point center() {
        return 0.5 * (vertices[0] + vertices[2]);
    }

    std::pair<Line, Line> diagonals() {
        return std::make_pair(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
    }
};

class Square : public Rectangle {
public:
    Square(const Point &p1, const Point &p2) : Rectangle(p1, p2, 1) {}
    virtual ~Square() {}
    
    Circle circumscribedCircle() {
        return Circle(this->center(), 0.5*dist(vertices[0], vertices[2]));    
    }

    Circle inscribedCircle() {
        return Circle(this->center(), 0.5/sqrt(2)*dist(vertices[0], vertices[2]));
    }
};

class Triangle : public Polygon {
public:
    Triangle (const std::vector<Point> &v) : Polygon(v) {}
    Triangle (const std::initializer_list<Point> &lst) : Polygon(lst) {}
    Triangle (const Point &a, const Point &b, const Point &c) : Polygon(a, b, c) {}
    virtual ~Triangle() {}
    
    Point centroid() {
        return (vertices[0]+vertices[1]+vertices[2])/3;
    }
    
    Point orthocenter() {
        Line l1(vertices[0], vertices[1]), l2(vertices[0], vertices[2]);
        Line l3(vertices[2], -1/l1.k), l4(vertices[1], -1/l2.k);
        return lineIntersection(l3, l4);
    }

    Point bisectorIntersection() {
        Line b1 = this->bisector(0), b2 = this->bisector(1);
        return lineIntersection(b1, b2);
    }

    Point perpendicularBisectorIntersection() {
        Line ab = Line(vertices[0], vertices[1]);
        Line bc = Line(vertices[1], vertices[2]);
        Line pb1 = Line(0.5*(vertices[0]+vertices[1]), -1/ab.k);
        Line pb2 = Line(0.5*(vertices[1]+vertices[2]), -1/bc.k);
        return lineIntersection(pb1, pb2);
    }

    Line EulerLine() {
        return Line(this->orthocenter(), this->centroid());
    }

    Line bisector(size_t i) {
        double a = dist(vertices[i], vertices[(i+1)%3]);
        double b = dist(vertices[i], vertices[(i+2)%3]);
        return Line(vertices[i], vertices[(i+1)%3]*b/(a+b) + vertices[(i+2)%3]*a/(a+b));
    }

    Circle circumscribedCircle() {
        Point center = perpendicularBisectorIntersection();
        return Circle(center, dist(vertices[0], center));
    }

    Circle inscribedCircle() {
        Point center = this->bisectorIntersection();
        return Circle(center, fabs(crossProduct(center-vertices[0], vertices[0]-vertices[1]))/dist(vertices[0], vertices[1]));
    }

    Circle ninePointsCircle() {
        Point orth = this->orthocenter();
        std::vector<Point> v(3);
        for (size_t i = 0; i < 3; i++) {
            v[i] = lineIntersection(Line(orth, vertices[i]), Line(vertices[(i+1)%3], vertices[(i+2)%3]));
        }
        Triangle orthotriangle(v);
        return orthotriangle.circumscribedCircle();
    }
};
