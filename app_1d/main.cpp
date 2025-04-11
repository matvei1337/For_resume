#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.h"

int main (int argc, char *argv[]) {
    QApplication app (argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *tool_bar = new QMenuBar (window);
    Window *graph_area = new Window (window);
    QAction *action;

    if (graph_area->parse_command_line (argc, argv)) {
        QMessageBox::warning (0, "Wrong input arguments!",
                            "Wrong input arguments!");
        return -1;
    }

    action = tool_bar->addAction ("&9 Previous func", graph_area, SLOT (previous_func ()));
    action->setShortcut (QString ("9"));

    action = tool_bar->addAction ("&0 Next func", graph_area, SLOT (next_func ()));
    action->setShortcut (QString ("0"));

    action = tool_bar->addAction ("&1 Change graph", graph_area, SLOT (change_graph ()));
    action->setShortcut (QString ("1"));

    action = tool_bar->addAction ("&2 Increase s", graph_area, SLOT (increase_s ()));
    action->setShortcut (QString ("2"));

    action = tool_bar->addAction ("&3 Decrease s", graph_area, SLOT (decrease_s ()));
    action->setShortcut (QString ("3"));

    action = tool_bar->addAction ("&4 Increase n", graph_area, SLOT (increase_n ()));
    action->setShortcut (QString ("4"));

    action = tool_bar->addAction ("&5 Decrease n", graph_area, SLOT (decrease_n ()));
    action->setShortcut (QString ("5"));

    action = tool_bar->addAction ("&6 Increase f", graph_area, SLOT (increase_f ()));
    action->setShortcut (QString("6"));

    action = tool_bar->addAction ("&7 Decrease f", graph_area, SLOT (decrease_f ()));
    action->setShortcut (QString("7"));

    action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
    action->setShortcut (QString ("Ctrl+X"));

    tool_bar->setMaximumHeight (30);

    window->setMenuBar (tool_bar);
    window->setCentralWidget (graph_area);
    window->setWindowTitle ("Graph");

    window->show ();
    app.exec ();
    delete window;
    return 0;
}
