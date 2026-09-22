#include <cmath>
#include <list>
#include <tuple>
#include <fstream>
#include <iostream>
#include <deque>
#include <memory>
#include <complex>
#include <functional>
#include <vector>
#include <algorithm>

#include "interp_internal.hh"

constexpr std::complex<double> I(0,1);

template <class T>
std::string to_string(T &d)
{
    std::ostringstream out;
    out << d;
    return out.str();
}

////////////////////////////////////////////////////////////////////////////////
class motion_base {
public:
    virtual ~motion_base() = default;
    virtual void straight_move(std::complex<double> end)=0;
    virtual void straight_rapid(std::complex<double> end)=0;
    virtual void circular_move(bool ccw,std::complex<double> center,
	std::complex<double> end)=0;
};

class motion_null:public motion_base {
public:
    void straight_move(std::complex<double> /*end*/) override  {}
    void straight_rapid(std::complex<double> /*end*/)  override {}
    void circular_move(bool /*ccw*/,std::complex<double> /*center*/,
	std::complex<double> /*end*/)  override {}
};

////////////////////////////////////////////////////////////////////////////////
class round_segment;
class straight_segment;

static constexpr double tolerance=1e-6;

class segment {
protected:
    std::complex<double> start, end;
    double finish{};
public:
    segment(double sz,double sx,double ez,double ex):start(sz,sx),end(ez,ex) {}
    segment(std::complex<double> s, std::complex<double> e):start(s),end(e) {}
    virtual ~segment() = default;
    typedef std::deque<double> intersections_t;
    virtual void intersection_z(double x, intersections_t &is)=0;
    virtual void crossings_up(double /*z*/, double /*xmin*/, int & /*n*/) {}
    virtual void draw(motion_base*)=0;
    virtual void offset(double)=0;
    virtual void intersect(segment*)=0;
    virtual void intersect_end(round_segment *p)=0;
    virtual void intersect_end(straight_segment *p)=0;
    virtual std::unique_ptr<segment> dup()=0;
    virtual double radius()=0;
    virtual bool monotonic() { return real(end-start)<=1e-3; }
    virtual void visit_x(std::function<void(double)> const &fn) { fn(imag(end)); }
    virtual void do_finish(segment * /*prev*/, segment * /*next*/) {}
    std::complex<double> &sp() { return start; }
    std::complex<double> &ep() { return end; }

    virtual void flip_imag() { start=conj(start); end=conj(end); }
    virtual void flip_real() { start=-conj(start); end=-conj(end); }
    virtual void rotate() {
	start=-start*I;
	end=-end*I;
    }
    virtual void move(std::complex<double> d) { start+=d; end+=d; }
    virtual void scale(double k) { start*=k; end*=k; finish*=k; }
    virtual double extent() {
	return std::max(std::max(std::abs(real(start)),std::abs(imag(start))),
	    std::max(std::abs(real(end)),std::abs(imag(end))));
    }
    friend class round_segment;

protected:
    std::complex<double> corner_finish(segment *prev, segment *next) {
	auto p=prev->dup();
	auto n=next->dup();

	p->offset(finish);
	n->offset(finish);
	if(real(p->end-n->start)>1e-2) {
	    finish=-finish;
	    p=prev->dup();
	    n=next->dup();
	    p->offset(finish);
	    n->offset(finish);
	}
	if(real(p->end-n->start)>1e-2)
	    throw(std::string("Corner finish failed at ") + to_string(prev->end));
	p->intersect(n.get());
	std::complex<double> center=(p->end+n->start)/2.0;
	p->offset(-finish);
	n->offset(-finish);
	start=prev->end=p->end;
	end=next->start=n->start;
	return center;
    }
};

class straight_segment:public segment {
public:
    straight_segment(double sx, double sz,double ex, double ez):
	segment(sz,sx,ez,ex) {}
    straight_segment(std::complex<double> s,std::complex<double> e):
	segment(s,e) {}
    void intersection_z(double x, intersections_t &is) override;
    void crossings_up(double z, double xmin, int &n) override;
    void draw(motion_base *out) override { out->straight_move(end); }
    void offset(double distance) override {
	std::complex<double> d=I*distance*(start-end)/std::abs(start-end);
	start+=d;
	end+=d;
    }
    void intersect(segment *p) override;
    void intersect_end(round_segment *p) override;
    void intersect_end(straight_segment *p) override;
    std::unique_ptr<segment> dup() override {
	return std::make_unique<straight_segment>(*this);
    }
    double radius() override { return std::abs(start-end); }
};

void straight_segment::intersection_z(double x, intersections_t &is)
{
    if(std::min(start.imag(),end.imag())>x+tolerance
	|| std::max(start.imag(),end.imag())<x-tolerance
    )
	return;
    if(std::abs(imag(start-end))<tolerance) {
	is.push_back(start.real());
	is.push_back(end.real());
    } else
	is.push_back((x-start.imag())
	    /(end.imag()-start.imag())*(end.real()-start.real())+start.real());
}

void straight_segment::crossings_up(double z, double xmin, int &n)
{
    double r0=real(start), r1=real(end);
    if(std::abs(r1-r0)<tolerance)
	return;
    double rmin=std::min(r0,r1), rmax=std::max(r0,r1);
    if(z<rmin-tolerance || z>rmax+tolerance)
	return;
    if(std::abs(z-rmax)<tolerance)
	return;
    double t=(z-r0)/(r1-r0);
    double im=imag(start)+t*(imag(end)-imag(start));
    if(im>xmin+tolerance)
	++n;
}

class round_segment:public segment {
protected:
    int ccw;
    std::complex<double> center;
public:
    round_segment(bool c, double sx, double sz,
	double cx, double cz, double ex, double ez):
	segment(sz,sx,ez,ex), ccw(c), center(cz,cx)
    {
    }
    round_segment(bool cc, std::complex<double> s,
	std::complex<double> c, std::complex<double> e):
	segment(s,e), ccw(cc), center(c)
    {
    }
    void intersection_z(double x,intersections_t &is) override;
    void crossings_up(double z, double xmin, int &n) override;
    void draw(motion_base *out) override { out->circular_move(ccw,center,end);}
    void offset(double distance) override {
	double factor=(std::abs(start-center)+(ccw? 1:-1)*distance)
	    /std::abs(start-center);
	start=factor*(start-center)+center;
	end=factor*(end-center)+center;
    }
    void intersect(segment *p) override;
    void intersect_end(round_segment *p) override;
    void intersect_end(straight_segment *p) override;
    std::complex<double> cp() { return center; }
    std::unique_ptr<segment> dup() override {
	return std::make_unique<round_segment>(*this);
    }
    double radius() override { return std::min(std::abs(start-end),std::abs(start-center)); }
    void flip_imag() override { ccw=!ccw; start=conj(start); end=conj(end);
	center=conj(center); }
    void flip_real() override { ccw=!ccw; start=-conj(start);
	end=-conj(end); center=-conj(center); }
    void rotate() override {
	start=-start*I;
	end=-end*I;
	center=-center*I;
    }
    bool monotonic() override {
	if(finish!=0)
	    return true;
	double entry=imag(start-center);
	double exit=imag(end-center);
	double dz=real(end-start);
	if(ccw)
	    return entry>=-1e-3 && exit>=-1e-3 && dz<=-1e-3;
	else
	    return entry<=1e-3 && exit<=1e-3 && dz<=-1e-3;
    }
    void visit_x(std::function<void(double)> const &fn) override;
    void move(std::complex<double> d) override { start+=d; center+=d; end+=d; }
    void scale(double k) override { start*=k; end*=k; center*=k; finish*=k; }
    double extent() override {
	return std::max(segment::extent(),
	    std::max(std::abs(real(center)),std::abs(imag(center))));
    }
private:
    bool on_segment(std::complex<double> p);
    friend class straight_segment;
};


inline bool round_segment::on_segment(std::complex<double> p)
{
    double radius=std::abs(start-center);
    if(!(radius>0))
	return true;
    // A complete circle: start and end coincide, any point on it counts.
    if(std::abs(start-end)<tolerance)
	return std::abs(std::abs(p-center)-radius)<=tolerance;
    /* Angular sweep, not two half planes: the intersection of the half
       planes is the wrong wedge for arcs spanning more than 180 degrees. */
    double two_pi=2*std::acos(-1.0);
    auto ang=[](std::complex<double> d) {
	return std::atan2(imag(d), real(d));
    };
    auto norm=[&](double a) {
	while(a<0)
	    a+=two_pi;
	while(a>=two_pi)
	    a-=two_pi;
	return a;
    };
    double a0=ang(start-center);
    double a1=ang(end-center);
    double ap=ang(p-center);
    double sweep=ccw? norm(a1-a0):norm(a0-a1);
    double d=ccw? norm(ap-a0):norm(a0-ap);
    // The angular tolerance scales with the radius, like the old test.
    return d<=sweep+tolerance/radius;
}

void round_segment::intersection_z(double x, intersections_t &is)
{
    // Decide on the distance from the centre, not on its square: an area
    // against a linear tolerance made a tangent scan line miss the arc
    double radius=std::abs(start-center);
    double dx=x-center.imag();
    if(std::abs(dx)>radius+tolerance)
	return;
    double s=(radius-dx)*(radius+dx);
    if(s<0)
	s=0;
    s=sqrt(s);
    std::complex<double> p1(real(center+s),x);
    if(on_segment(p1))
	is.push_back(real(p1));
    std::complex<double> p2(real(center-s),x);
    if(on_segment(p2))
	is.push_back(real(p2));
}

void round_segment::crossings_up(double z, double xmin, int &n)
{
    double radius=std::abs(start-center);
    double dz=z-real(center);
    if(std::abs(dz)>radius+tolerance)
	return;
    double s=(radius-dz)*(radius+dz);
    if(s<0)
	s=0;
    s=sqrt(s);
    std::complex<double> p1(z, imag(center)+s);
    std::complex<double> p2(z, imag(center)-s);
    if(imag(p1)>xmin+tolerance && on_segment(p1))
	++n;
    if(imag(p2)>xmin+tolerance && on_segment(p2))
	++n;
}

void round_segment::visit_x(std::function<void(double)> const &fn)
{
    double r=std::abs(start-center);
    if (!(r>0)) {
	fn(imag(end));
	return;
    }
    auto ang=[](std::complex<double> p, std::complex<double> c) {
	auto d=p-c;
	return std::atan2(imag(d), real(d));
    };
    double a0=ang(start, center);
    double a1=ang(end, center);
    double pi=std::acos(-1.0);
    double two_pi=2*pi;
    double sweep=ccw? (a1-a0):(a0-a1);
    if (sweep<=0)
	sweep+=two_pi;
    double ext[2]={pi/2, -pi/2};
    double dist[2], xs[2];
    int n=0;
    for (double ae : ext) {
	double d=ccw? (ae-a0):(a0-ae);
	while (d<0)
	    d+=two_pi;
	while (d>=two_pi)
	    d-=two_pi;
	if (d>1e-9 && d<sweep-1e-9) {
	    dist[n]=d;
	    xs[n]=imag(center)+r*std::sin(ae);
	    n++;
	}
    }
    if (n==2 && dist[1]<dist[0]) {
	std::swap(dist[0], dist[1]);
	std::swap(xs[0], xs[1]);
    }
    for (int i=0; i<n; i++)
	fn(xs[i]);
    fn(imag(end));
}

////////////////////////////////////////////////////////////////////////////////
void straight_segment::intersect(segment *p)
{
    p->intersect_end(this);
}

void round_segment::intersect(segment *p)
{
    p->intersect_end(this);
}

void straight_segment::intersect_end(straight_segment *p)
{
    // correct end of p and start of this
    auto len=std::abs(start-end);
    if(!(len>tolerance)) {
	start=p->end;
	return;
    }
    auto rot=conj(start-end)/len;
    auto ps=(p->start-end)*rot;
    auto pe=(p->end-end)*rot;
    if(std::abs(imag(ps-pe))<tolerance) {
	/* Parallel: do not throw (a closing shoulder often hits this).
	   Close the gap at the midpoint so offset/D can continue. */
	auto mid=(p->end+start)/2.0;
	start=p->end=mid;
	return;
    }
    auto f=imag(ps)/imag(ps-pe);
    auto is=(ps+f*(pe-ps))/rot+end;
    /* A near-parallel pair can push the vertex arbitrarily far and look
       like a bulge.  Keep the join near the original ends. */
    if(std::abs(is-start)>10*len && std::abs(is-p->end)>10*std::abs(p->end-p->start)) {
	auto mid=(p->end+start)/2.0;
	start=p->end=mid;
	return;
    }
    start=p->end=is;
}

void straight_segment::intersect_end(round_segment *p)
{
    // correct end of p and start of this
    // (arc followed by a straight
    if(std::abs(start-p->end)<tolerance)
	return;

    auto rot=conj(start-end)/std::abs(start-end);
    auto pe=(p->end-end)*rot;
    auto pc=(p->center-end)*rot;

    double b=norm(pc-pe)-imag(pc)*imag(pc);
    if(b<0) {
	b=0;
    } else
	b=sqrt(b);

    auto s1=(real(pc)+b)/rot+end;
    auto s2=(real(pc)-b)/rot+end;

    if(std::abs(start-s1)<std::abs(start-s2))
	start=p->end=s1;
    else
	start=p->end=s2;
}

void round_segment::intersect_end(straight_segment *p)
{
    // correct end of p and start of this
    // (straight followed by an arc)
    if(std::abs(start-p->end)<tolerance)
	return;

    auto rot=conj(p->start-p->end)/std::abs(p->start-p->end);
    auto pe=(end-p->end)*rot;
    auto pc=(center-p->end)*rot;

    double b=norm(pc-pe)-imag(pc)*imag(pc);
    if(b<0) {
	b=0;
    } else
	b=sqrt(b);

    auto s1=(real(pc)+b)/rot+p->end;
    auto s2=(real(pc)-b)/rot+p->end;

    if(std::abs(start-s1)<std::abs(start-s2))
	start=p->end=s1;
    else
	start=p->end=s2;
}

void round_segment::intersect_end(round_segment *p)
{
    // correct end of p and start of this
    auto a=std::abs(start-center);
    auto b=std::abs(p->start-p->center);
    auto c=std::abs(center-p->center);
    auto cosB=(c*c+a*a-b*b)/2.0/a/c;
    if(std::abs(cosB)>1) {
	// circles do not intersect: trim to the closest points on the
	// line of centers so each end stays on its own circle
	auto u=(p->center-center)/c;
	p->end=p->center-b*u;
	start=center+a*u;
	return;
    }
    double cosB2=cosB*cosB;
    std::complex<double> rot(cosB,sqrt(1-cosB2));
    auto is=rot*(p->center-center)/c*a+center;
    p->end=start=is;
}

////////////////////////////////////////////////////////////////////////////////
// Corner finishes

class fillet_segment:public round_segment {
public:
    fillet_segment(double d,std::complex<double> l):round_segment(0,l,l,l) {
	finish=d;
    }

    void do_finish(segment *prev, segment *next) override {
	center=corner_finish(prev,next);
	ccw=imag(start-center)>0 || imag(end-center)>0;
	finish=0;
    }
};

class chamfer_segment:public straight_segment {
public:
    chamfer_segment(double d,std::complex<double> l):straight_segment(l,l) {
	finish=d;
    }

    void do_finish(segment *prev, segment *next) override {
	corner_finish(prev,next);
    }
};


////////////////////////////////////////////////////////////////////////////////

class scaled_motion:public motion_base {
    motion_base *orig;
    double k;
    std::complex<double> location;
    bool located{false};
public:
    scaled_motion(motion_base *motion,double scale):orig(motion),k(scale) {}
    void straight_move(std::complex<double> end) override {
	orig->straight_move(k*end);
	location=end; located=true;
    }
    void straight_rapid(std::complex<double> end) override {
	orig->straight_rapid(k*end);
	location=end; located=true;
    }
    void circular_move(bool ccw,std::complex<double> center,
	std::complex<double> end) override
    {
	// Drop an arc going nowhere, judged on the normalised profile
	if(located && std::abs(end-location)<tolerance)
	    return;
	orig->circular_move(ccw,k*center,k*end);
	location=end; located=true;
    }
};

template <int swap>
class swapped_motion:public motion_base {
    motion_base *orig;
public:
    swapped_motion(motion_base *motion):orig(motion) {
    }
    void straight_move(std::complex<double> end) override {
	switch(swap) {
	case 0: orig->straight_move(end); break;
	case 1: orig->straight_move(-conj(end)); break;
	case 2: orig->straight_move(conj(end)); break;
	case 3: orig->straight_move(-end); break;
	case 4: orig->straight_move(I*end); break;
	case 5: orig->straight_move(conj(I*end)); break;
	case 6: orig->straight_move(-conj(I*end)); break;
	case 7: orig->straight_move(-I*end); break;
	}
    }
    void straight_rapid(std::complex<double> end) override {
	switch(swap) {
	case 0: orig->straight_rapid(end); break;
	case 1: orig->straight_rapid(-conj(end)); break;
	case 2: orig->straight_rapid(conj(end)); break;
	case 3: orig->straight_rapid(-end); break;
	case 4: orig->straight_rapid(I*end); break;
	case 5: orig->straight_rapid(conj(I*end)); break;
	case 6: orig->straight_rapid(-conj(I*end)); break;
	case 7: orig->straight_rapid(-I*end); break;
	}
    }
    void circular_move(bool ccw,std::complex<double> center,
	std::complex<double> end) override
    {
	switch(swap) {
	case 0: orig->circular_move(ccw,center,end); break;
	case 1: orig->circular_move(!ccw,-conj(center),-conj(end)); break;
	case 2: orig->circular_move(!ccw,conj(center),conj(end)); break;
	case 3: orig->circular_move(ccw,-center,-end); break;
	case 4: orig->circular_move(ccw,I*center,I*end); break;
	case 5: orig->circular_move(!ccw,conj(I*center),conj(I*end)); break;
	case 6: orig->circular_move(!ccw,-conj(I*center),-conj(I*end)); break;
	case 7: orig->circular_move(ccw,-I*center,-I*end); break;
	}
    }
};

////////////////////////////////////////////////////////////////////////////////
class g7x:public  std::list<std::unique_ptr<segment>> {
    double delta{0.5};
    std::complex<double> escape{0.3, 0.3};
    int flip_state{};
    double unit{1};
    std::unique_ptr<motion_base> unscaled;
private:
    void scan_rough(int cycle, motion_base *out);
    void add_distance(double distance);
    bool in_stock(double z, double x);

    /* Rotate profile by 90 degrees */
    void rotate() {
	for(auto &p : *this)
	    p->rotate();
	flip_state^=4;
    }

    /* Change the direction of the profile to have Z and X decreasing */
    void swap() {
	double dir_x=imag(front()->ep()-back()->ep());
	double dir_z=real(front()->ep()-back()->ep());
	if(dir_x>0) {
	    for(auto &p : *this)
		p->flip_imag();
	    flip_state^=2;
	}
	if(dir_z<0) {
	    for(auto &p : *this)
		p->flip_real();
	    flip_state^=1;
	}
    }

    /* The thresholds here are plain numbers, so the profile needs a known
       size for them to mean anything.  Scale it into a canonical range, by a
       power of two so the scaling itself is exact, and undo that on output.
    */
    void normalise() {
	double size=0;
	for(auto &p : *this)
	    size=std::max(size,p->extent());
	if(!(size>0) || !std::isfinite(size))
	    return;
	int exponent;
	std::frexp(size,&exponent);
	unit=std::ldexp(1.0,6-exponent);	// scaled size lands in [32,64)
	if(unit==1)
	    return;
	for(auto &p : *this)
	    p->scale(unit);
    }

    std::unique_ptr<motion_base> motion(motion_base *out) {
	unscaled=std::make_unique<scaled_motion>(out,1/unit);
	out=unscaled.get();
	switch(flip_state) {
	case 0: return std::make_unique<swapped_motion<0>>(out);
	case 1: return std::make_unique<swapped_motion<1>>(out);
	case 2: return std::make_unique<swapped_motion<2>>(out);
	case 3: return std::make_unique<swapped_motion<3>>(out);
	case 4: return std::make_unique<swapped_motion<4>>(out);
	case 5: return std::make_unique<swapped_motion<5>>(out);
	case 6: return std::make_unique<swapped_motion<6>>(out);
	case 7: return std::make_unique<swapped_motion<7>>(out);
	}
	throw(std::string("This can't happen"));
    }

    void monotonic(const char *axis="Z") {
	if(real(front()->ep()-front()->sp())>0) {
	    front()->sp().real(real(front()->ep()));
	}
	for(auto &p : *this) {
	    if(!p->monotonic())
		throw(std::string("Not monotonic in ") + axis);
	}
    }

    void do_finish() {
	for(auto h=++begin(); h!=--end(); ) {
	    auto p(h); --p;
	    (*h)->do_finish((*p).get(),(*(++h)).get());
	}
	for(auto p=begin(); p!=end();)
	    if((*p)->radius()<1e-3)
		p=erase(p);
	    else
		++p;
    }

    bool should_rotate_paths() {
	g7x paths(*this);
	paths.pop_front();
	paths.swap();
	for(auto &path : paths) {
	    if(!path->monotonic())
		return true;
	}
	return false;
    }

public:
    g7x() = default;
    g7x(g7x const &other)
          : std::list<std::unique_ptr<segment>>(),
            delta(other.delta),
            escape(other.escape),
            flip_state(other.flip_state),
            unit(other.unit) {
	for(const auto &p : other)
	    emplace_back(p->dup());
    }

    /*
	x,z	from where the distance is computed, also a rapid to this
		location between each pass.
	d	Starting distance
	e	Ending distance
	p	Number of passes to go from d to e
    */
    void do_g70(motion_base *out, double x, double z, double d, double e,
	int p, CUTTER_COMP *cutter_comp_side, bool rapid_entry=true
    )
    {
	front()->sp()=std::complex<double>(z,x);
	normalise();
	d*=unit;
	e*=unit;

	bool rotated=should_rotate_paths();
	if(rotated)
	    rotate();
	swap();
	do_finish();
	monotonic(rotated ? "X" : "Z");
	auto swapped_out=motion(out);

	for(int pass=p; pass>0; pass--) {
	    double distance=(pass-1)*(d-e)/p+e;
	    g7x paths(*this);
	    paths.add_distance(distance);

	    auto comp=*cutter_comp_side;
	    *cutter_comp_side=CUTTER_COMP::OFF;
	    swapped_out->straight_rapid(paths.front()->sp());
	    *cutter_comp_side=comp;
	    if (rapid_entry)
		swapped_out->straight_rapid(paths.front()->ep());
	    else
		paths.front()->draw(swapped_out.get());
	    paths.pop_front();
	    for(const auto &path : paths)
		path->draw(swapped_out.get());
	}
    }

    /*
	cycle	0 Complete cut including pockets.
		1 Do not cut any pockets.
		2 Only cut after first pocket (pick up where 1 stopped)
	x,z	From where the cutting is done.
	d	Final distance to profile
	i	Increment of cutting (depth of cut)
	r	Distance to retract
    */
    void do_g71(motion_base *out, int cycle, double x, double z,
	double u, double w,
	double d, double i, double r, bool do_rotate=false,
	bool do_envelope=false
    ) {
	front()->sp()=std::complex<double>(z,x);
	normalise();
	u*=unit;
	w*=unit;
	d*=unit;
	i*=unit;
	r*=unit;
	if(do_rotate)
	    rotate();
	swap();
	do_finish();
	monotonic(do_rotate ? "X" : "Z");
	add_distance(d);
	std::complex<double> displacement(w,u);
	for(auto &p : *this)
	    p->move(displacement);
	auto swapped_out=motion(out);

	delta=std::max(i,2*tolerance);
	escape=r*std::complex<double>{1,1};

	swapped_out->straight_rapid(front()->sp());
	// Close the profile only if it really is open, not a few ulp short
	if(imag(back()->ep())+tolerance<imag(front()->sp())) {
	    auto ep=back()->ep();
	    ep.imag(imag(front()->sp()));
	    emplace_back(std::make_unique<straight_segment>(
		back()->ep(),ep
	    ));
	}
	scan_rough(cycle, swapped_out.get());
	if(do_envelope) {
	    /* Follow the offset profile once so a native G71/G72 does not
	       leave steps of one depth of cut on tapers and arcs.  The path
	       is already offset by D, so this leaves a uniform finishing
	       allowance.  Cycles that skip pockets must not do this. */
	    swapped_out->straight_rapid(front()->sp());
	    for(auto &p : *this)
		p->draw(swapped_out.get());
	}
    }

    void do_g72(motion_base *out, int cycle, double x, double z,
	double u, double w,
	double d, double i, double r, bool do_envelope=false
    ) {
	do_g71(out, cycle, x, z, u, w, d, i, r, true, do_envelope);
    }
};


bool g7x::in_stock(double z, double x)
{
    int n=0;
    for(auto &p : *this)
	p->crossings_up(z, x, n);
    return (n%2)==0;
}

void g7x::scan_rough(int cycle, motion_base *out)
{
    if(empty())
	return;
    double x_top=imag(front()->sp());
    double z0=real(front()->sp());
    double x_min=x_top, x_max=x_top;
    double z_lo=z0, z_hi=z0;
    auto track_x=[&](double xv) {
	x_min=std::min(x_min, xv);
	x_max=std::max(x_max, xv);
    };
    for(auto &p : *this) {
	/* visit_x also sees the inside of arcs, so the deepest pass reaches
	   the real bottom of a radius, not just its end points. */
	p->visit_x(track_x);
	z_lo=std::min(z_lo, std::min(real(p->sp()), real(p->ep())));
	z_hi=std::max(z_hi, std::max(real(p->sp()), real(p->ep())));
    }
    if(!(delta>0) || x_min>=x_top-tolerance)
	return;
    /* Rapids travel at x_clear, above every point of the profile, so no
       rapid ever moves in Z at a cutting X. */
    double x_clear=std::max(x_top, x_max);

    auto unique_eps=[](std::vector<double> &v, double eps) {
	std::sort(v.begin(), v.end());
	std::vector<double> o;
	for(double t : v) {
	    if(o.empty() || std::abs(t-o.back())>eps)
		o.push_back(t);
	}
	v.swap(o);
    };

    struct Interval { double za, zb; };
    auto intervals_at=[&](double x)->std::vector<Interval> {
	std::vector<double> zs;
	for(auto &p : *this) {
	    segment::intersections_t is;
	    p->intersection_z(x, is);
	    for(double t : is)
		zs.push_back(t);
	}
	zs.push_back(z0);
	zs.push_back(z_lo);
	zs.push_back(z_hi);
	unique_eps(zs, tolerance);
	std::vector<Interval> cuts;
	for(size_t i=0; i+1<zs.size(); ++i) {
	    double a=zs[i], b=zs[i+1];
	    if(std::abs(b-a)<2*tolerance)
		continue;
	    if(in_stock(0.5*(a+b), x))
		cuts.push_back({a, b});
	}
	auto dist=[&](const Interval &iv) {
	    return std::min(std::abs(iv.za-z0), std::abs(iv.zb-z0));
	};
	std::sort(cuts.begin(), cuts.end(),
	    [&](const Interval &A, const Interval &B) { return dist(A)<dist(B); });
	/* The interval holding the start face is the outer envelope that
	   G7x.1 already roughed; pockets are everything else. */
	auto envelope=[&](const Interval &iv) {
	    return iv.za-tolerance<=z0 && z0<=iv.zb+tolerance;
	};
	if(cycle==1) {
	    cuts.erase(std::remove_if(cuts.begin(), cuts.end(),
		[&](const Interval &iv) { return !envelope(iv); }), cuts.end());
	} else if(cycle==2) {
	    cuts.erase(std::remove_if(cuts.begin(), cuts.end(), envelope),
		cuts.end());
	}
	return cuts;
    };

    double passes=(x_top-x_min)/delta+8;
    if(!(passes<100000))
	passes=100000;
    int max_pass=(int)passes;
    if(max_pass<2)
	max_pass=2;

    for(int pass=1; pass<=max_pass; ++pass) {
	double x=x_top-pass*delta;
	bool last=false;
	if(x<=x_min+tolerance) {
	    x=x_min;
	    last=true;
	}
	auto cuts=intervals_at(x);
	for(auto &iv : cuts) {
	    double z_enter, z_exit;
	    if(std::abs(iv.za-z0)<=std::abs(iv.zb-z0)) {
		z_enter=iv.za;
		z_exit=iv.zb;
	    } else {
		z_enter=iv.zb;
		z_exit=iv.za;
	    }
	    /* Travel at the clearance X and plunge at the entry Z: no rapid
	       ever moves to a cutting X while crossing the part. */
	    out->straight_rapid(std::complex<double>(z_enter, x_clear));
	    out->straight_rapid(std::complex<double>(z_enter, x));
	    out->straight_move(std::complex<double>(z_exit, x));
	    double distance=std::abs(z_exit-z_enter);
	    double er=std::max(std::abs(real(escape)), tolerance);
	    double esc=std::min(1.0, distance/2/er);
	    std::complex<double> up=std::complex<double>(z_exit, x)+esc*escape;
	    out->straight_rapid(up);
	    out->straight_rapid(std::complex<double>(real(up), x_clear));
	}
	out->straight_rapid(std::complex<double>(z0, x_clear));
	if(last)
	    break;
    }
}

void g7x::add_distance(double distance) {
    auto v1=front()->ep()-front()->sp();
    auto v2=back()->sp()-front()->sp();
    auto angle=v1/v2;
    if(imag(angle)<0)
	distance=-distance;
    auto of(std::move(front()));
    pop_front();
    if(empty())
	throw(std::string("Profile disappeared while applying stock allowance"));
    /* Offset in steps of at most half the shortest segment.  Track what is
       left, not what has been added, so the last step lands exactly instead
       of leaving a residue worth another pass.  A step that no longer
       shortens the remainder ends the loop rather than spinning.
    */
    double remaining=std::abs(distance);
    double sign=distance<0? -1.0:1.0;
    while(remaining>tolerance) {
	if(empty())
	    throw(std::string("Profile disappeared while applying stock allowance"));
	double step=1e9;
	for(auto &p : *this)
	    step=std::min(step,p->radius()/2);
	if(!(step>0))
	    break;
	if(step>=remaining-tolerance)
	    step=remaining;
	remaining-=step;
	double max_distance=sign*step;

	for(auto p=begin(); p!=end();) {
	    (*p)->offset(max_distance);
	    if((*p)->radius()<1e-3)
		p=erase(p);
	    else
		++p;
	}

	if(size()<2)
	    break;

	for(auto p=begin(); p!=--end(); ++p) {
	    auto n=p; ++n;
	    if(real((*p)->ep()-(*n)->sp())>1e-2) {
		/* Convex-corner gap: a fillet of radius ~|step| is correct.
		   A much larger arc is the #2939 bulge; join with a straight. */
		auto s((*p)->dup());
		auto e((*n)->dup());
		s->offset(-max_distance);
		e->offset(-max_distance);
		auto center=(s->ep()+e->sp())/2.0;
		double arcr=std::abs((*p)->ep()-center);
		if(arcr>std::abs(max_distance)*2+1e-3) {
		    emplace(n, std::make_unique<straight_segment>(
			(*p)->ep(),(*n)->sp()));
		} else {
		    emplace(n, std::make_unique<round_segment>(
			distance>0,(*p)->ep(),center,(*n)->sp()));
		}
		++p;
	    }
	}

	if(size()<2)
	    break;

	for(auto p=begin(); p!=--end(); ++p) {
	    if(!(*p)->monotonic()) {
		auto pp=p; --pp;
		if((*p)->radius()<(*pp)->radius())
		    erase(p);
		else
		    erase(pp);
		p=begin();
	    }
	    auto n=p; ++n;
	    if((*p)->radius()<1e-3) {
		erase(p);
		p=begin();
	    } else
		(*p)->intersect(n->get());
	}
    }
    for(auto p=begin(); p!=--end(); ++p) {
	auto n=p; ++n;
	auto gap=std::abs((*p)->ep()-(*n)->sp());
	if(gap>1e-3)
	    throw(std::string("Profile junction gap too large: ")
		+to_string((*p)->ep()));
	auto mid=((*p)->ep()+(*n)->sp())/2.0;
	(*p)->ep()=(*n)->sp()=mid;
    }

    of->ep()=front()->sp();
    push_front(std::move(of));
}

#ifndef IGNORE_LINUXCNC
////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <cstdio>
#include <string>
#include "rs274ngc.hh"
#include "rs274ngc_interp.hh"
#include "interp_internal.hh"
#include "interp_queue.hh"

/* Cycle words on a G7x block must not leak into convert_straight/arc
   (P as arc turns, R as radius-format, U/W as extra-axis motion). */
static void g7x_clear_cycle_words(block_pointer block)
{
    block->a_flag=false;
    block->b_flag=false;
    block->c_flag=false;
    block->d_flag=false;
    block->e_flag=false;
    block->i_flag=false;
    block->j_flag=false;
    block->k_flag=false;
    block->p_flag=false;
    block->p_number=-1;
    block->q_flag=false;
    block->q_number=-1;
    block->r_flag=false;
    block->u_flag=false;
    block->v_flag=false;
    block->w_flag=false;
}

class motion_machine:public motion_base {
    Interp *interp;
    setup_pointer settings;
    block_pointer block;
    void prepare_xz(std::complex<double> end) {
	g7x_clear_cycle_words(block);
	block->x_flag=1;
	block->x_number=imag(end);
	block->z_flag=1;
	block->z_number=real(end);
    }
public:
    motion_machine(Interp *i, setup_pointer s, block_pointer b):
	interp(i), settings(s), block(b) { }

    void straight_move(std::complex<double> end) override {
	prepare_xz(end);
	int r=interp->convert_straight(G_1, block, settings);
	if(r!=INTERP_OK)
	    throw(r);
    }

    void straight_rapid(std::complex<double> end) override  {
	prepare_xz(end);
	int r=interp->convert_straight(G_0, block, settings);
	if(r!=INTERP_OK)
	    throw(r);
    }
    void circular_move(bool ccw,std::complex<double> center,
	std::complex<double> end
    ) override {
	prepare_xz(end);
	block->i_flag=1;
	block->i_number=imag(center);
	block->k_flag=1;
	block->k_number=real(center);

	// scaled_motion has already dropped the empty ones
	int r=interp->convert_arc(ccw? G_3:G_2, block, settings);
	if(r!=INTERP_OK)
	    throw(r);
    }
};

class switch_settings {
    Interp *interp;
    setup_pointer settings;
    DISTANCE_MODE saved_distance_mode, saved_ijk_distance_mode;
    read_function_pointer read_a, read_c, read_u, read_w;
public:
    switch_settings(Interp *i,setup_pointer s):interp(i), settings(s)
    {
	saved_distance_mode=settings->distance_mode;
	settings->distance_mode=DISTANCE_MODE::ABSOLUTE;
	saved_ijk_distance_mode=settings->ijk_distance_mode;
	settings->ijk_distance_mode=DISTANCE_MODE::ABSOLUTE;
	read_a=interp->_readers[(int)'a'];
	read_c=interp->_readers[(int)'c'];
	read_u=interp->_readers[(int)'u'];
	read_w=interp->_readers[(int)'w'];
	interp->_readers[(int)'a']=Interp::default_readers[(int)'a'];
	interp->_readers[(int)'c']=Interp::default_readers[(int)'c'];
	interp->_readers[(int)'u']=Interp::default_readers[(int)'u'];
	interp->_readers[(int)'w']=Interp::default_readers[(int)'w'];
    }
    ~switch_settings() {
	settings->distance_mode=saved_distance_mode;
	settings->ijk_distance_mode=saved_ijk_distance_mode;
	interp->_readers[(int)'a']=read_a;
	interp->_readers[(int)'c']=read_c;
	interp->_readers[(int)'u']=read_u;
	interp->_readers[(int)'w']=read_w;
    }
    DISTANCE_MODE ijk_distance_mode() { return saved_ijk_distance_mode; }
    DISTANCE_MODE distance_mode() { return saved_distance_mode; }
};

int Interp::convert_g7x(int mode,
      block_pointer block,     //!< pointer to a block of RS274 instructions
      setup_pointer settings)  //!< pointer to machine settings
{
    int cycle=mode;
    int subcycle=cycle%10;
    cycle/=10;
    bool fanuc = (mode == G_71_3 || mode == G_70_3);

    /* Fanuc parameter line: G71.3 U(Δd) R(e) stores depth/retract only. */
    if (mode == G_71_3 && !block->p_flag && !block->q_flag) {
	if (block->u_flag) {
	    settings->g71_3_delta = block->u_number;
	    settings->g71_3_have_delta = true;
	}
	if (block->r_flag) {
	    settings->g71_3_retract = block->r_number;
	    settings->g71_3_have_retract = true;
	}
	return INTERP_OK;
    }

    if (fanuc) {
	CHKS((!block->p_flag || !block->q_flag),
	    _("G%d.%d execute line requires P and Q sequence numbers"),
	    cycle, subcycle);
    } else if (!block->q_flag) {
	ERS("G7x.x  requires a Q word");
    }

    if(settings->cutter_comp_side != CUTTER_COMP::OFF && cycle!=70)
	ERS("G%d.%d cannot be used with cutter compensation enabled",
	    cycle, subcycle);
    if(settings->plane!=CANON_PLANE::XZ)
	ERS("G%d.%d can only be used in XZ plane (G18)",
	    cycle, subcycle);

    switch_settings old(this,settings);

    auto original_block=*block;

    double x=settings->current_x;
    double z=settings->current_z;
    double start_x=x;
    double start_z=z;
    if(old.distance_mode()==DISTANCE_MODE::INCREMENTAL) {
	if(block->x_flag)
	    x+=block->x_number;
	if(block->z_flag)
	    z+=block->z_number;
    } else {
	if(block->x_flag)
	    x=block->x_number;
	if(block->z_flag)
	    z=block->z_number;
    }
    original_block.x_number=x;
    original_block.z_number=z;
    g7x_clear_cycle_words(block);

    /* On ERS/CHP, put the interpreter pose back so a failed cycle does not
       leave current_x/z / motion_mode at a half-built profile point. */
    struct RestorePose {
	setup_pointer s;
	double cx, cz;
	int motion;
	bool armed;
	RestorePose(setup_pointer st):
	    s(st), cx(st->current_x), cz(st->current_z),
	    motion(st->motion_mode), armed(true) {}
	void dismiss() { armed=false; }
	~RestorePose() {
	    if (armed) {
		s->current_x=cx;
		s->current_z=cz;
		s->motion_mode=motion;
	    }
	}
    } pose(settings);

    if (!fanuc || original_block.x_flag || original_block.z_flag) {
	auto cutter_comp=settings->cutter_comp_side;
	settings->cutter_comp_side=CUTTER_COMP::OFF;
	int error=convert_straight(G_0, block, settings);
	settings->cutter_comp_side=cutter_comp;
	if(error!=INTERP_OK)
	    return error;
    }

    g7x path;
    std::complex<double> start(z,x);
    bool type_ii=false;
    bool first_profile_block=true;
    int first_motion_mode=-1;
    bool have_profile_f=false;
    double profile_f=0;

    auto apply_params = [&]() -> int {
	for(int n=0; n<settings->parameter_occurrence; n++)
	    settings->parameters[settings->parameter_numbers[n]]=
		settings->parameter_values[n];

	for(int n=0; n<settings->named_parameter_occurrence; n++)
	    CHP(store_named_param(&_setup, settings->named_parameters[n],
		settings->named_parameter_values[n]
	    ));
	settings->named_parameter_occurrence = 0;
	return INTERP_OK;
    };

    auto append_block = [&]() -> int {
	if(block->g_modes[GM_MOTION]!=-1)
	    settings->motion_mode=block->g_modes[GM_MOTION];
	else if (settings->motion_mode > 30) {
	    if (fanuc)
		ERS(_("G%d.%d first profile block must contain G0, G1, G2 or G3"),
		    cycle, subcycle);
	    else
		settings->motion_mode = 10;
	}

	if (first_profile_block &&
	    (block->x_flag || block->z_flag || block->u_flag || block->w_flag)) {
	    type_ii = block->z_flag || block->w_flag;
	    first_motion_mode = settings->motion_mode;
	    first_profile_block = false;
	}

	if (fanuc && cycle == 70 && block->f_flag && !have_profile_f) {
	    have_profile_f = true;
	    profile_f = block->f_number;
	}

	std::complex<double> end(start);
	std::complex<double> center(0,0);

	if(old.distance_mode()==DISTANCE_MODE::INCREMENTAL)
	    end=0;

	if(block->u_flag) {
	    if(old.distance_mode()==DISTANCE_MODE::INCREMENTAL)
		ERS("G7x error: Cannot use U in incremental mode (G91)");
	    if(block->x_flag)
		ERS("G7x error: Cannot use U and X in the same block");
	    auto uu=block->u_number;
	    if(settings->lathe_diameter_mode)
		uu/=2;
	    end.imag(end.imag()+uu);
	} else if(block->x_flag)
	    end.imag(block->x_number);

	if(block->w_flag) {
	    if(old.distance_mode()==DISTANCE_MODE::INCREMENTAL)
		ERS("G7x error: Cannot use W in incremental mode (G91)");
	    if(block->z_flag)
		ERS("G7x error: Cannot use W and Z in the same block");
	    end.real(end.real()+block->w_number);
	} else if(block->z_flag)
	    end.real(block->z_number);

	if(block->i_flag) center.imag(block->i_number);
	if(block->k_flag) center.real(block->k_number);

	if(old.distance_mode()==DISTANCE_MODE::INCREMENTAL)
	    end+=start;

	if(start!=end) {
	    switch(settings->motion_mode) {
	    case 0:
	    case 10:
		path.emplace_back(std::make_unique<straight_segment>(
		    start, end
		));
		break;
	    case 20:
	    case 30:
		if(block->r_flag) {
		    if(block->i_flag || block->k_flag)
			ERS("G7X error: both R and I or K flag used for arc");
		    double rr=block->r_number;
		    center=(start+end)/2.0;
		    CHKS((rr*rr+tolerance<norm(end-start)/4.0),
			_("G7X error: arc radius too small to reach the end point"));
		    auto dd=I*sqrt((rr*rr-norm(end-start)/4)/norm(end-start))
			*(end-start);
		    if(settings->motion_mode==30)
			center+=dd;
		    else
			center-=dd;
		} else {
		    if(!block->i_flag && !block->k_flag)
			ERS("G7X error: either I or K must be present for arc");
		    if(old.ijk_distance_mode()==DISTANCE_MODE::INCREMENTAL)
			center+=start;
		}

		double tpx;
		double tpy;
		double tpz;
		double tAA_p;
		double tBB_p;
		double tCC_p;
		double tu_p;
		double tv_p;
		double tw_p;
		CHP(find_ends(block, settings, &tpx, &tpy, &tpz, &tAA_p, &tBB_p, &tCC_p, &tu_p, &tv_p, &tw_p));

		if(!block->r_flag){
		    double center1, center2;
		    int turn;
		    double radius_tolerance = (settings->length_units == CANON_UNITS_INCHES) ? RADIUS_TOLERANCE_INCH : RADIUS_TOLERANCE_MM;
		    double spiral_abs_tolerance = (settings->length_units == CANON_UNITS_INCHES) ? settings->center_arc_radius_tolerance_inch : settings->center_arc_radius_tolerance_mm;
		    CHP(arc_data_ijk((settings->motion_mode==30)? G_3 : G_2, settings->plane, settings->current_z, settings->current_x, tpz, tpx,
				     (old.ijk_distance_mode() == DISTANCE_MODE::ABSOLUTE),
				     (block->k_flag)? block->k_number : 0.0, (block->i_flag)? block->i_number : 0.0, block->p_flag? round_to_int(block->p_number) : 1,
				     &center1, &center2, &turn, radius_tolerance, spiral_abs_tolerance, SPIRAL_RELATIVE_TOLERANCE));
		}

		path.emplace_back(std::make_unique<round_segment>(
		    settings->motion_mode==30, start, center, end
		));
		break;
	    }
	    if(settings->motion_mode==0 || settings->motion_mode==10
		|| settings->motion_mode==20 || settings->motion_mode==30
	    ) {
		if(block->a_flag && block->c_flag)
		    ERS("G7X error: Both A and C parameters on a corner");
		if(block->a_flag)
		    path.emplace_back(std::make_unique<fillet_segment>(
			block->a_number, end
		    ));
		if(block->c_flag)
		    path.emplace_back(std::make_unique<chamfer_segment>(
			block->c_number, end
		    ));
	    }

	    settings->current_x=imag(end);
	    settings->current_z=real(end);
	    start=end;
	}
	return INTERP_OK;
    };

    if (fanuc) {
	CHKS((settings->file_pointer == NULL),
	    _("G%d.%d P/Q profile requires an open program file (not MDI)"),
	    cycle, subcycle);
	int n_start = round_to_int(original_block.p_number);
	int n_end = round_to_int(original_block.q_number);
	CHKS((n_start < 0 || n_end < 0),
	    _("G%d.%d P and Q must be non-negative sequence numbers"),
	    cycle, subcycle);

	struct RestoreFile {
	    FILE *fp;
	    long pos;
	    int *seqn;
	    int seq;
	    RestoreFile(setup_pointer s):
		fp(s->file_pointer),
		pos(s->file_pointer ? ftell(s->file_pointer) : -1),
		seqn(&s->sequence_number),
		seq(s->sequence_number) {}
	    ~RestoreFile() {
		if (fp && pos >= 0)
		    fseek(fp, pos, SEEK_SET);
		if (seqn)
		    *seqn = seq;
	    }
	} restore(settings);

	/* Speculatively apply # assignments while reading the profile, then
	   put parameters back. Lines between G71.3 and N(P) still run once
	   after the cycle; lines inside P-Q are skipped. */
	struct RestoreParams {
	    setup_pointer s;
	    std::vector<double> numbered;
	    parameter_map named0;
	    parameter_map named_local;
	    int level;
	    RestoreParams(setup_pointer st):
		s(st),
		numbered(st->parameters,
		    st->parameters + interp_param_global::RS274NGC_MAX_PARAMETERS),
		named0(st->sub_context[0].named_params),
		named_local(st->sub_context[st->call_level].named_params),
		level(st->call_level) {}
	    void restore() {
		std::copy(numbered.begin(), numbered.end(), s->parameters);
		s->sub_context[0].named_params = named0;
		s->sub_context[level].named_params = named_local;
	    }
	    ~RestoreParams() { restore(); }
	} restore_params(settings);

	/* Search forward from this cycle first: a duplicate N before G71.3
	   must not be taken as the profile.  If the range is not found there
	   (G70.3 is normally placed after the profile), rewind and search
	   from the start of the program, as Fanuc does. */
	CHKS((restore.pos < 0),
	    _("G%d.%d could not read the program to find N%d"),
	    cycle, subcycle, n_start);

	bool found_p=false;
	bool found_q=false;
	int scan_motion_mode=settings->motion_mode;
	long p_line_pos=-1, q_line_pos=-1;
	/* Also record where the range starts and ends, so a later G70.3 with
	   the same P/Q numbers can read exactly what G71.3 used. */
	auto do_scan=[&](long start_pos) -> int {
	    found_p=false;
	    found_q=false;
	    p_line_pos=-1;
	    q_line_pos=-1;
	    CHKS((fseek(settings->file_pointer, start_pos, SEEK_SET) != 0),
		_("G%d.%d could not read the program to find N%d"),
		cycle, subcycle, n_start);
	    bool collecting=false;
	    char raw_line[LINELEN];
	    char line[LINELEN];
	    for(;;) {
		long line_pos=ftell(settings->file_pointer);
		if (fgets(raw_line, LINELEN, settings->file_pointer) == NULL)
		    break;
		settings->parameter_occurrence = 0;
		settings->named_parameter_occurrence = 0;
		if (strlen(raw_line) == (LINELEN - 1))
		    ERS("G7X error: line too long while reading P-Q profile");
		for (int index = (int)strlen(raw_line) - 1;
		     index >= 0 && isspace(static_cast<unsigned char>(raw_line[index]));
		     index--)
		    raw_line[index] = 0;
		rs274ngc_strlcpy(line, raw_line, LINELEN);
		CHP(close_and_downcase(line));
		/* Block delete, like read_text: only drop an N hidden on a
		   '/' line when block delete is actually enabled. */
		if (line[0] == 0 ||
		    (line[0] == '/' && GET_BLOCK_DELETE()) ||
		    (line[0] == '%' && line[1] == 0))
		    continue;

		CHP(parse_line(line, block, settings));
		/* Apply # assignments on every scanned line, including those
		   between G71.3 and N(P), so the profile sees the new values. */
		CHP(apply_params());

		if (!collecting && block->n_number == n_start) {
		    collecting = true;
		    found_p = true;
		    p_line_pos = line_pos;
		}
		if (collecting) {
		    CHP(append_block());
		    if (block->n_number == n_end) {
			found_q = true;
			q_line_pos = line_pos;
			break;
		    }
		}
	    }
	    return INTERP_OK;
	};

	auto reset_scan=[&]() {
	    /* Undo the speculative # of a failed pass before collecting
	       again. */
	    restore_params.restore();
	    path.clear();
	    start=std::complex<double>(z,x);
	    settings->current_x=imag(start);
	    settings->current_z=real(start);
	    settings->motion_mode=scan_motion_mode;
	    first_profile_block=true;
	    type_ii=false;
	    first_motion_mode=-1;
	    have_profile_f=false;
	    profile_f=0;
	};

	/* G70.3 normally sits after the profile that G71.3 already found, so
	   reuse that exact range when the numbers match.  Otherwise search:
	   G70.3 from the start of the program (Fanuc), G71.3 forward from
	   the cycle; fall back to the other order if the range is missing. */
	bool cached = (cycle==70) &&
	    settings->g7x_profile_valid &&
	    (strcmp(settings->g7x_profile_file, settings->filename)==0) &&
	    settings->g7x_profile_p==n_start &&
	    settings->g7x_profile_q==n_end;
	bool done=false;
	if (cached) {
	    CHP(do_scan(settings->g7x_profile_p_pos));
	    done = found_p && found_q;
	    if (!done)
		reset_scan();
	}
	if (!done) {
	    CHP(do_scan(cycle==70 ? 0 : restore.pos));
	    if (!found_p || !found_q) {
		reset_scan();
		CHP(do_scan(cycle==70 ? restore.pos : 0));
	    }
	}
	if (cycle==71 && found_p && found_q && p_line_pos>=0 && q_line_pos>=0) {
	    rs274ngc_strlcpy(settings->g7x_profile_file, settings->filename, LINELEN);
	    settings->g7x_profile_p_pos=p_line_pos;
	    settings->g7x_profile_q_pos=q_line_pos;
	    settings->g7x_profile_p=n_start;
	    settings->g7x_profile_q=n_end;
	    settings->g7x_profile_valid=true;
	}
	CHKS(!found_p, _("G%d.%d P sequence number N%d not found in the program"),
	    cycle, subcycle, n_start);
	CHKS(!found_q, _("G%d.%d Q sequence number N%d not found after N%d"),
	    cycle, subcycle, n_end, n_start);
	*block = original_block;
	g7x_clear_cycle_words(block);
    } else {
	auto exit_call_level=settings->call_level;
	/* Q was cleared from `block` before the optional start rapid so it
	   cannot leak into convert_straight; the subroutine number lives
	   on original_block. */
	CHP(read((std::string("O")+std::to_string(static_cast<int>(original_block.q_number))+" CALL").c_str()));
	for(;;) {
	    if(block->o_name!=NULL)
		CHP(convert_control_functions(block, settings));
	    if(settings->call_level==exit_call_level)
		break;
	    CHP(apply_params());
	    CHP(append_block());
	    CHP(read());
	}
    }
    if(path.size()<=1) {
	if (fanuc)
	    ERS(_("G%d.%d profile between P and Q is empty or too short"),
		cycle, subcycle);
	pose.dismiss();
	return INTERP_OK;
    }

    double d=0, e=0, i=1, r=0.5, u=0, w=0;
    int p = 1;
    if(original_block.d_flag) d=original_block.d_number_float;
    if(original_block.e_flag) e=original_block.e_number;
    if(original_block.i_flag) i=original_block.i_number;
    if(original_block.p_flag) p=static_cast<int>(original_block.p_number);
    if(original_block.r_flag) r=original_block.r_number;
    if(original_block.u_flag) u=original_block.u_number;
    if(original_block.w_flag) w=original_block.w_number;
    if(original_block.x_flag) x=original_block.x_number;
    if(original_block.z_flag) z=original_block.z_number;

    if (fanuc && cycle == 71) {
	if (original_block.i_flag)
	    i = original_block.i_number;
	else if (settings->g71_3_have_delta)
	    i = settings->g71_3_delta;
	else
	    ERS(_("G71.3 requires depth of cut (U on the parameter line or I on the execute line)"));
	if (original_block.r_flag)
	    r = original_block.r_number;
	else if (settings->g71_3_have_retract)
	    r = settings->g71_3_retract;
	else
	    r = (settings->length_units == CANON_UNITS_INCHES) ? 0.020 : 0.5;
	/* Finish allowance U is in the same units as X: diameter in G7, radius in G8. */
	u = original_block.u_flag ? original_block.u_number : 0.0;
	if (original_block.u_flag && settings->lathe_diameter_mode)
	    u /= 2.0;
	w = original_block.w_flag ? original_block.w_number : 0.0;
	d = 0;
    }
    if (fanuc && cycle == 70) {
	d = 0;
	e = 0;
	p = 1;
	u = 0;
	w = 0;
	if (!original_block.f_flag && have_profile_f) {
	    settings->feed_rate = profile_f;
	    enqueue_SET_FEED_RATE(profile_f);
	}
    }

    settings->current_x=start_x;
    settings->current_z=start_z;

    if(i<=0)
	ERS("G7X error: I must be greater than zero.");

    auto dfront=path.front()->ep()-path.front()->sp();
    auto dback=path.back()->ep()-path.back()->sp();
    motion_machine motion(this, settings, block);
    try {
	switch(cycle) {
	case 70:
	    path.do_g70(&motion,x,z,d,e,p,&settings->cutter_comp_side,
		!(fanuc && (first_motion_mode == 10 || first_motion_mode == 20 ||
			    first_motion_mode == 30)));
	    break;
	case 71:
	    if (fanuc && !type_ii) {
		auto it=path.begin();
		if (it != path.end())
		    ++it; /* skip A to A'; ignore a later closing segment */
		double xcur=0;
		int dir=0;
		bool have_x=false;
		if (it != path.end()) {
		    xcur=imag((*it)->sp());
		    have_x=true;
		}
		for (; it != path.end(); ++it) {
		    bool failed=false;
		    (*it)->visit_x([&](double xv) {
			if (failed || !have_x)
			    return;
			if (std::abs(xv-xcur)<tolerance) {
			    xcur=xv;
			    return;
			}
			int s=xv>xcur? 1:-1;
			if (!dir)
			    dir=s;
			else if (s!=dir)
			    failed=true;
			xcur=xv;
		    });
		    if (failed)
			ERS(_("G71.3 Type I profile is not monotonic in X"));
		}
	    }
	    if(std::abs(real(dback))>0 && imag(dfront)*(x-imag(start))<0) {
		std::complex<double> end{real(start),x};
		path.emplace_back(std::make_unique<straight_segment>(
		    start, end
		));
	    }
	    path.do_g71(&motion, fanuc ? (type_ii ? 0 : 1) : subcycle,
		x,z,u,w,d,i,r, false, !fanuc && subcycle==0);
	    break;
	case 72:
	    if(std::abs(imag(dback))>0 && real(dfront)*(z-real(start))<0) {
		std::complex<double> end{z,imag(start)};
		path.emplace_back(std::make_unique<straight_segment>(
		    start, end
		));
	    }
	    path.do_g72(&motion,subcycle,x,z,u,w,d,i,r,
		!fanuc && subcycle==0);
	    break;
	}
    } catch(std::string &s) {
	ERS("G7X error: %s", s.c_str());
	return INTERP_ERROR;
    } catch(int err) {
	return err;
    }

    if (fanuc) {
	auto cutter_comp=settings->cutter_comp_side;
	settings->cutter_comp_side=CUTTER_COMP::OFF;
	g7x_clear_cycle_words(block);
	block->x_flag=1;
	block->x_number=start_x;
	block->z_flag=1;
	block->z_number=start_z;
	int error=convert_straight(G_0, block, settings);
	settings->cutter_comp_side=cutter_comp;
	if(error!=INTERP_OK)
	    return error;
	settings->g7x_skip_n_start = round_to_int(original_block.p_number);
	settings->g7x_skip_n_end = round_to_int(original_block.q_number);
	settings->g7x_skip_active = false;
    }

    pose.dismiss();
    return INTERP_OK;
}
#endif
