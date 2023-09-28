QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    src/addbead.cpp \
    src/addbeadtype.cpp \
    src/addbrush.cpp \
    src/addshell.cpp \
    src/main.cpp \
    src/mainwindow.cpp \
    src/tipswindow.cpp

HEADERS += \
    src/addbead.h \
    src/addbeadtype.h \
    src/addbrush.h \
    src/addshell.h \
    src/mainwindow.h \
    src/tipswindow.h

FORMS += \
    src/addbead.ui \
    src/addbeadtype.ui \
    src/addbrush.ui \
    src/addshell.ui \
    src/mainwindow.ui \
    src/tipswindow.ui

OBJECTS_DIR = build
MOC_DIR = build
RCC_DIR = build
UI_DIR = build


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
