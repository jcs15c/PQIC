class Unit : public F_params {
public:
    double eval( Point x ) { return 1.0; };
};

class Cx : public F_params {
public:
    double eval( Point x ) { return x[0]; };
};

class Cy : public F_params {
public:
    double eval( Point x ) { return x[1]; };
};

class Mx2 : public F_params {
public:
    double eval( Point x ) { return x[0]*x[0]; };
};

class My2 : public F_params {
public:
    double eval( Point x ) { return x[1]*x[1]; };
};

class Mxy : public F_params {
public:
    double eval( Point x ) { return x[0]*x[1]; };
};

////////////////////////////////////
class Gx4 : public F_params {
public:
    Gx4() { };
    Gx4( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0]*x[0]*x[0] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx2y2 : public F_params {
public:
    Gx2y2() { };
    Gx2y2( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0]*x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx3y : public F_params {
public:
    Gx3y() { };
    Gx3y( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0]*x[0]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx3 : public F_params {
public:
    Gx3() { };
    Gx3( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0]*x[0] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx2y : public F_params {
public:
    Gx2y() { };
    Gx2y( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gy4 : public F_params {
public:
    Gy4() { };
    Gy4( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[1]*x[1]*x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gxy3 : public F_params {
public:
    Gxy3() { };
    Gxy3( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[1]*x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gxy2 : public F_params {
public:
    Gxy2() { };
    Gxy2( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx2 : public F_params {
public:
    Gx2() { };
    Gx2( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[0] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gy3 : public F_params {
public:
    Gy3() { };
    Gy3( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[1]*x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gxy : public F_params {
public:
    Gxy() { };
    Gxy( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gy2 : public F_params {
public:
    Gy2() { };
    Gy2( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[1]*x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gy : public F_params {
public:
    Gy() { };
    Gy( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[1] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

class Gx : public F_params {
public:
    Gx() { };
    Gx( Phi* phi00 ) { phi0 = phi00; };
    double eval( Point x ) { 
        return x[0] / gradient( x, phi0 );
    }; 
    Phi* phi0;
};

