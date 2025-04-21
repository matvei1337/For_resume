#include <QPainter>
#include <stdio.h>
#include "window.h"

#define L2G(X,Y) (l2g((X), (Y), min_y, max_y))

static double normalize(double value, double abs_min, double abs_max){
    if(value - abs_min < 0) return 0;
    if(abs_max - value < 0) return 1;
    if(std::fabs(abs_max - abs_min) < EPS) return 1;
    return std::fabs((value - abs_min) / (abs_max - abs_min));
}

static void rgb_for_f(int& R, int& G, int& B, double k){
    R = k * 128;
    G = k * 0;
    B = k * 128;
}

static void rgb_for_pf(int& R, int& G, int& B, double k){
    R = k * 0;
    G = k * 255;
    B = k * 255;
}

static void rgb_for_residual(int& R, int& G, int& B, double k){
    R = k * 255;
    G = k * 255;
    B = k * 0;
}



Window::Window(QWidget *parent): QWidget(parent){
    widget = parent;
}

Window::~Window(){
    if(args) delete[] args;
    if(tid) delete[] tid;
}



QSize Window::minimumSizeHint() const{
    return QSize(100, 100);
}

QSize Window::sizeHint() const{
    return QSize(1000, 1000);
}

QPointF Window::l2g(double x_loc, double y_loc, double y_min, double y_max){
    double x_gl = (x_loc - data.a) / (data.b - data.a) * width();
    double y_gl = (y_max - y_loc) / (y_max - y_min) * height();
    return QPointF(x_gl, y_gl);
}



int Window::parse_command_line(int argc, char* argv[]){
    if(argc != 13) return -1;

    data.read_data(argv);
    data.realloc_data();

    args = new Args[data.p];
    tid = new pthread_t[data.p];
    program_name = argv[0];
    func_id = data.m - 1;

    for(int i = 0; i < data.p; i++){
        args[i].copy_data(data);
        args[i].k = i;
        args[i].mutex = &mutex;
        args[i].cond = &cond;
    }

    next_func();
    return 0;
}

void Window::select_f(){
    double (*f)(double, double) = data.f;
    switch(func_id){
        case 0:
            f_name = "f(x, y) = 1";
            f = f0;
            break;
        case 1:
            f_name = "f(x, y) = x";
            f = f1;
            break;
        case 2:
            f_name = "f(x, y) = y";
            f = f2;
            break;
        case 3:
            f_name = "f(x, y) = x + y";
            f = f3;
            break;
        case 4:
            f_name = "f(x, y) = sqrt(x*x + y*y)";
            f = f4;
            break;
        case 5:
            f_name = "f(x, y) = x*x + y*y";
            f = f5;
            break;
        case 6:
            f_name = "f(x, y) = exp(x*x - y*y)";
            f = f6;
            break;
        case 7:
            f_name = "f(x, y) = 1 / (25(x*x + y*y) + 1)";
            f = f7;
            break;
    }
    data.f = f;
}

void Window::triangle(QPainter* painter, Point p1, Point p2, Point p3, QColor color){
    QPainterPath path;
    double min_y = data.c, max_y = data.d;

    path.moveTo(L2G(p1.x, p1.y));
    path.lineTo(L2G(p2.x, p2.y));
    path.lineTo(L2G(p3.x, p3.y));
    path.lineTo(L2G(p1.x, p1.y));

    painter->setPen(Qt::NoPen);
    painter->fillPath(path, QBrush(color));
}

void Window::draw_f(QPainter* painter){
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my; double (*f)(double, double) = data.f;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, k_rgb;
    int R, G, B;
    Point p1, p2, p3;
    double abs_min = f(a + hx/3.0, c + (2.0/3.0) * hy);
    double abs_max = abs_min;
    data.find_min_max(abs_min, abs_max);
    for(int i = 0; i < mx; ++i){
        for (int j = 0; j < my; ++j){
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            value = f(a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy);
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_f(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));

            p2 = {a + (i + 1)*hx, c + j*hy};
            value = f(a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy);
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_f(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_str[1000];
    double f_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_str, "max_abs(f) = %.2e", f_max);
    painter->setPen("white");
    painter->drawText(5, 40, max_str);
    draw_txt(painter);
}

void Window::draw_Pf(QPainter* painter){
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, k_rgb;
    int R, G, B;
    Point p1, p2, p3;
    double abs_min = data.pf({a + hx/3.0, c + (2.0/3.0) * hy});
    double abs_max = abs_min;

    data.pfind_min_max(abs_min, abs_max);
    for(int i = 0; i < mx; ++i){
        for(int j = 0; j < my; ++j){
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            value = data.pf({a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy});
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_pf(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));

            p2 = {a + (i + 1)*hx, c + j*hy};
            value = data.pf({a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy});
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_pf(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_str[1000];
    double pf_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_str, "max_abs(Pf) = %.2e", pf_max);
    painter->setPen("blue");
    painter->drawText(5, 40, max_str);
    draw_txt(painter);
}


void Window::draw_residual(QPainter *painter){
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my; double (*f)(double, double) = data.f;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, f_val, pf_val, k_rgb;
    int R, G, B;
    Point p1, p2, p3;
    double abs_min = std::fabs(f(a + hx/3.0, c + (2.0/3.0) * hy) - data.pf({a + hx/3.0, c + (2.0/3.0) * hy}));
    double abs_max = abs_min;
    data.residual_min_max(abs_min, abs_max);

    for(int i = 0; i < mx; ++i){
        for(int j = 0; j < my; ++j){
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            f_val = f(a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy);
            pf_val = data.pf({a + (i + 1.0/3.0)*hx, c + (j + 2.0/3.0) * hy});

            value = std::fabs(f_val - pf_val);
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_residual(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));

            p2 = {a + (i + 1)*hx, c + j*hy};
            f_val = f(a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy);
            pf_val = data.pf({a + (i + 2.0/3.0)*hx, c + (j + 1.0/3.0) * hy});

            value = std::fabs(f_val - pf_val);
            k_rgb = normalize(value, abs_min, abs_max);
            rgb_for_residual(R, G, B, k_rgb);
            triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_str[1000];
    double res_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_str, "max_abs(Pf - f) = %.2e", res_max);
    painter->setPen("yellow");
    painter->drawText(5, 40, max_str);
    draw_txt(painter);
}

void Window::draw_txt(QPainter* painter){
    painter->drawText(5, 20, f_name);
    std::string txt1 = "a = " + std::to_string(data.a);
    std::string txt2 = "b = " + std::to_string(data.b);
    std::string txt3 = "c = " + std::to_string(data.c);
    std::string txt4 = "d = " + std::to_string(data.d);
    std::string txt5 = "nx = " + std::to_string(data.nx);
    std::string txt6 = "ny = " + std::to_string(data.ny);
    std::string txt7 = "mx = " + std::to_string(data.mx);
    std::string txt8 = "my = " + std::to_string(data.my);
    std::string txt9 = "p = " + std::to_string(p_count);;
    painter->drawText(5, 60, txt1.c_str());
    painter->drawText(5, 80, txt2.c_str());
    painter->drawText(5, 100, txt3.c_str());
    painter->drawText(5, 120, txt4.c_str());
    painter->drawText(5, 140, txt5.c_str());
    painter->drawText(5, 160, txt6.c_str());
    painter->drawText(5, 180, txt7.c_str());
    painter->drawText(5, 200, txt8.c_str());
    painter->drawText(5, 220, txt9.c_str());
}

bool Window::msr_ready(){
    for(int i = 0; i < data.p; i++){
        if(!args[i].ready) return false;
    }
    return true;
}

void Window::msr_wait(){
    if(msr_ready()){
        printf("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n\
                It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",program_name, task, args[0].res1, args[0].res2, args[0].res3, args[0].res4,
                args[0].t1, args[0].t2, args[0].its, data.eps, func_id, data.nx, data.ny, data.p);
        update();
    }
    else{
        QTimer::singleShot(200, this, &Window::msr_wait);
    }
}




void Window::next_func(){
    if(msr_ready()){
        func_id = (func_id + 1) % 8;
        select_f();

        double abs_min, abs_max;
        f_max = data.find_min_max(abs_min, abs_max);

        data.realloc_data();

        for(int i = 0; i < data.p; i++){
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        if(threads_created) pthread_cond_broadcast(&cond);
        else{
            threads_created = true;
            for(int i = 0; i < data.p; i++){
                pthread_create(&args[i].tid, 0, solution, args + i);
            }
        }
        QTimer::singleShot(200, this, &Window::msr_wait);
    }
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}

void Window::next_graph(){
    n_graph = (n_graph + 1) % 3;
    update();
}

void Window::inc_n(){
    if(msr_ready()){
        data.nx *= 2;
        data.ny *= 2;
        data.realloc_data();

        for(int i = 0; i < data.p; i++){
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }
        pthread_cond_broadcast(&cond);
        QTimer::singleShot(200, this, &Window::msr_wait);
    }
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}

void Window::dec_n(){
    if(msr_ready()){
        data.nx /= 2;
        data.ny /= 2;
        data.realloc_data();

        for(int i = 0; i < data.p; i++){
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond);
        QTimer::singleShot(200, this, &Window::msr_wait);
    } 
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}

void Window::inc_p(){
    if(msr_ready()){
        p_count += 1;
        data.realloc_data();

        for(int i = 0; i < data.p; i++){
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond);
        QTimer::singleShot(200, this, &Window::msr_wait);
    }
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}

void Window::dec_p(){
    if(msr_ready()){
        p_count -= 1;
        data.realloc_data();

        for(int i = 0; i < data.p; i++){
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond);
        QTimer::singleShot(200, this, &Window::msr_wait);
    }
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}

void Window::inc_s(){
    double len_x = data.b - data.a;
    double len_y = data.d - data.c;

    data.a += len_x / 4;
    data.b -= len_x / 4;
    data.c += len_y / 4;
    data.d -= len_y / 4;

    update();
}

void Window::dec_s(){
    double len_x = data.b - data.a;
    double len_y = data.d - data.c;

    data.a -= len_x / 2;
    data.b += len_x / 2;
    data.c -= len_y / 2;
    data.d += len_y / 2;

    update();
}

void Window::inc_m(){
    data.mx *= 2;
    data.my *= 2;

    double abs_min, abs_max;
    f_max = data.find_min_max(abs_min, abs_max);

    update();
}

void Window::dec_m(){
    data.mx /= 2;
    data.my /= 2;

    double abs_min, abs_max;
    f_max = data.find_min_max(abs_min, abs_max);

    update();
}

void Window::close(){
    if(msr_ready()) widget->close();
    else QMessageBox::warning(0, "Waiting msr!!!", "Calculation not complete!");
}



void Window::paintEvent(QPaintEvent*){
    QPainter painter(this);
    if(n_graph == 0){
        draw_f(&painter);
    }
    else if(n_graph == 1){
        draw_Pf(&painter);
    }
    else{
        draw_residual(&painter);
    }
}
