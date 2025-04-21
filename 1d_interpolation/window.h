#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget {
    Q_OBJECT

private:
    int func_id;
    const char* f_name;
    double a;
    double b;
    int n;
    int k;
    int n_graph = 0;
    int s = 0;
    int p = 0;
    double (*f) (double);
    double (*df) (double);

public:
    Window (QWidget *parent);

    QSize minimumSizeHint () const;
    QSize sizeHint () const;

    int parse_command_line (int argc, char *argv[]);
    void select_f();

    void draw_Pf(QPainter* painter, const QVector<double>& x, const QVector<double>& y);
    void draw_Sf(QPainter* painter, const QVector<double>& x, const QVector<double>& y);
    void draw_f(QPainter* painter);
    void draw_residual_Pf(QPainter* painter, const QVector<double>& x, const QVector<double>& y);
    void draw_residual_Sf(QPainter* painter, const QVector<double>& x, const QVector<double>& y);

     QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
public slots:
    void previous_func();
    void next_func();
    void change_graph();
    void increase_s();
    void decrease_s();
    void increase_n();
    void decrease_n();
    void increase_f();
    void decrease_f();

protected:
    void paintEvent (QPaintEvent *event);
};

#endif
