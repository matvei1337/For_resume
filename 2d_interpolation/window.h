#pragma once
#include <QtWidgets/QtWidgets>
#include "solver.h"
#include "gui.h"

class Window : public QWidget{
    Q_OBJECT
private:
    QWidget* widget;
    int func_id;
    const char* f_name;
    const char* program_name;
    int n_graph = 0;
    double f_max;
    int p_count = 0;

    Args* args = nullptr;
    pthread_t* tid = nullptr;
    Gui_data data;
    const int task = 3;
    bool threads_created = false;
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

public:
    Window(QWidget* parent);
    ~Window();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    QPointF l2g(double x_loc, double y_loc, double y_min, double y_max);

    int parse_command_line(int argc, char *argv[]);
    void select_f();
    void triangle(QPainter* painter, Point p1, Point p2, Point p3, QColor color);
    void draw_f(QPainter* painter);
    void draw_Pf(QPainter* painter);
    void draw_residual(QPainter* painter);
    void draw_txt(QPainter* painter);
    bool msr_ready();
    void msr_wait();

public slots:
    void next_func();
    void next_graph();
    void inc_s();
    void dec_s();
    void inc_n();
    void dec_n();
    void inc_p();
    void dec_p();
    void inc_m();
    void dec_m();
    void close();

protected:
    void paintEvent(QPaintEvent *event);
};