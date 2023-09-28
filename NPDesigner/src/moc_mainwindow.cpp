/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[56];
    char stringdata0[854];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 28), // "on_newBeadTypeButton_clicked"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 24), // "on_newBeadButton_clicked"
QT_MOC_LITERAL(4, 66, 25), // "on_newShellButton_clicked"
QT_MOC_LITERAL(5, 92, 25), // "on_newBrushButton_clicked"
QT_MOC_LITERAL(6, 118, 16), // "recieveNewNPBead"
QT_MOC_LITERAL(7, 135, 1), // "x"
QT_MOC_LITERAL(8, 137, 1), // "y"
QT_MOC_LITERAL(9, 139, 1), // "z"
QT_MOC_LITERAL(10, 141, 6), // "beadID"
QT_MOC_LITERAL(11, 148, 8), // "doUpdate"
QT_MOC_LITERAL(12, 157, 20), // "recieveNewNPBeadType"
QT_MOC_LITERAL(13, 178, 11), // "std::string"
QT_MOC_LITERAL(14, 190, 13), // "hamakerFileIn"
QT_MOC_LITERAL(15, 204, 12), // "surfaceDirIn"
QT_MOC_LITERAL(16, 217, 8), // "radiusIn"
QT_MOC_LITERAL(17, 226, 18), // "surfacePotentialIn"
QT_MOC_LITERAL(18, 245, 12), // "surfFactorIn"
QT_MOC_LITERAL(19, 258, 12), // "coreFactorIn"
QT_MOC_LITERAL(20, 271, 10), // "ljCutoffIn"
QT_MOC_LITERAL(21, 282, 20), // "correctionOverrideIn"
QT_MOC_LITERAL(22, 303, 15), // "recieveNewShell"
QT_MOC_LITERAL(23, 319, 11), // "innerRadius"
QT_MOC_LITERAL(24, 331, 11), // "outerRadius"
QT_MOC_LITERAL(25, 343, 8), // "ljCutoff"
QT_MOC_LITERAL(26, 352, 15), // "recieveNewBrush"
QT_MOC_LITERAL(27, 368, 14), // "brushOccupancy"
QT_MOC_LITERAL(28, 383, 15), // "brushRadialDist"
QT_MOC_LITERAL(29, 399, 10), // "beadTypeID"
QT_MOC_LITERAL(30, 410, 11), // "forceAttach"
QT_MOC_LITERAL(31, 422, 19), // "updateBeadTypeTable"
QT_MOC_LITERAL(32, 442, 15), // "updateBeadTable"
QT_MOC_LITERAL(33, 458, 7), // "updateM"
QT_MOC_LITERAL(34, 466, 1), // "m"
QT_MOC_LITERAL(35, 468, 9), // "numPoints"
QT_MOC_LITERAL(36, 478, 11), // "updateTheta"
QT_MOC_LITERAL(37, 490, 5), // "theta"
QT_MOC_LITERAL(38, 496, 1), // "j"
QT_MOC_LITERAL(39, 498, 23), // "on_actionQuit_triggered"
QT_MOC_LITERAL(40, 522, 20), // "on_findUADir_clicked"
QT_MOC_LITERAL(41, 543, 23), // "on_saveNPButton_clicked"
QT_MOC_LITERAL(42, 567, 23), // "on_actionTips_triggered"
QT_MOC_LITERAL(43, 591, 20), // "updateGraphicsWindow"
QT_MOC_LITERAL(44, 612, 23), // "on_updateTables_clicked"
QT_MOC_LITERAL(45, 636, 28), // "on_beadTypeTable_cellChanged"
QT_MOC_LITERAL(46, 665, 3), // "row"
QT_MOC_LITERAL(47, 669, 6), // "column"
QT_MOC_LITERAL(48, 676, 21), // "on_recenterNP_clicked"
QT_MOC_LITERAL(49, 698, 18), // "updateBindingRadii"
QT_MOC_LITERAL(50, 717, 26), // "on_autoBounds_stateChanged"
QT_MOC_LITERAL(51, 744, 4), // "arg1"
QT_MOC_LITERAL(52, 749, 28), // "on_lowerBoundLine_textEdited"
QT_MOC_LITERAL(53, 778, 28), // "on_upperBoundLine_textEdited"
QT_MOC_LITERAL(54, 807, 22), // "on_actionNew_triggered"
QT_MOC_LITERAL(55, 830, 23) // "on_actionLoad_triggered"

    },
    "MainWindow\0on_newBeadTypeButton_clicked\0"
    "\0on_newBeadButton_clicked\0"
    "on_newShellButton_clicked\0"
    "on_newBrushButton_clicked\0recieveNewNPBead\0"
    "x\0y\0z\0beadID\0doUpdate\0recieveNewNPBeadType\0"
    "std::string\0hamakerFileIn\0surfaceDirIn\0"
    "radiusIn\0surfacePotentialIn\0surfFactorIn\0"
    "coreFactorIn\0ljCutoffIn\0correctionOverrideIn\0"
    "recieveNewShell\0innerRadius\0outerRadius\0"
    "ljCutoff\0recieveNewBrush\0brushOccupancy\0"
    "brushRadialDist\0beadTypeID\0forceAttach\0"
    "updateBeadTypeTable\0updateBeadTable\0"
    "updateM\0m\0numPoints\0updateTheta\0theta\0"
    "j\0on_actionQuit_triggered\0"
    "on_findUADir_clicked\0on_saveNPButton_clicked\0"
    "on_actionTips_triggered\0updateGraphicsWindow\0"
    "on_updateTables_clicked\0"
    "on_beadTypeTable_cellChanged\0row\0"
    "column\0on_recenterNP_clicked\0"
    "updateBindingRadii\0on_autoBounds_stateChanged\0"
    "arg1\0on_lowerBoundLine_textEdited\0"
    "on_upperBoundLine_textEdited\0"
    "on_actionNew_triggered\0on_actionLoad_triggered"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      26,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  144,    2, 0x08 /* Private */,
       3,    0,  145,    2, 0x08 /* Private */,
       4,    0,  146,    2, 0x08 /* Private */,
       5,    0,  147,    2, 0x08 /* Private */,
       6,    5,  148,    2, 0x08 /* Private */,
      12,    8,  159,    2, 0x08 /* Private */,
      22,    6,  176,    2, 0x08 /* Private */,
      26,    4,  189,    2, 0x08 /* Private */,
      31,    0,  198,    2, 0x08 /* Private */,
      32,    0,  199,    2, 0x08 /* Private */,
      33,    2,  200,    2, 0x08 /* Private */,
      36,    3,  205,    2, 0x08 /* Private */,
      39,    0,  212,    2, 0x08 /* Private */,
      40,    0,  213,    2, 0x08 /* Private */,
      41,    0,  214,    2, 0x08 /* Private */,
      42,    0,  215,    2, 0x08 /* Private */,
      43,    0,  216,    2, 0x08 /* Private */,
      44,    0,  217,    2, 0x08 /* Private */,
      45,    2,  218,    2, 0x08 /* Private */,
      48,    0,  223,    2, 0x08 /* Private */,
      49,    0,  224,    2, 0x08 /* Private */,
      50,    1,  225,    2, 0x08 /* Private */,
      52,    1,  228,    2, 0x08 /* Private */,
      53,    1,  231,    2, 0x08 /* Private */,
      54,    0,  234,    2, 0x08 /* Private */,
      55,    0,  235,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Int, QMetaType::Bool,    7,    8,    9,   10,   11,
    QMetaType::Void, 0x80000000 | 13, 0x80000000 | 13, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Int,   14,   15,   16,   17,   18,   19,   20,   21,
    QMetaType::Void, 0x80000000 | 13, 0x80000000 | 13, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Double,   14,   15,   23,   24,   17,   25,
    QMetaType::Void, QMetaType::Double, QMetaType::Double, QMetaType::Int, QMetaType::Bool,   27,   28,   29,   30,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Double, QMetaType::Double, QMetaType::Int,   34,   35,
    QMetaType::Double, QMetaType::Double, QMetaType::Double, QMetaType::Double,   37,   34,   38,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   46,   47,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   51,
    QMetaType::Void, QMetaType::QString,   51,
    QMetaType::Void, QMetaType::QString,   51,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MainWindow *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->on_newBeadTypeButton_clicked(); break;
        case 1: _t->on_newBeadButton_clicked(); break;
        case 2: _t->on_newShellButton_clicked(); break;
        case 3: _t->on_newBrushButton_clicked(); break;
        case 4: _t->recieveNewNPBead((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4])),(*reinterpret_cast< bool(*)>(_a[5]))); break;
        case 5: _t->recieveNewNPBeadType((*reinterpret_cast< std::string(*)>(_a[1])),(*reinterpret_cast< std::string(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3])),(*reinterpret_cast< float(*)>(_a[4])),(*reinterpret_cast< float(*)>(_a[5])),(*reinterpret_cast< float(*)>(_a[6])),(*reinterpret_cast< float(*)>(_a[7])),(*reinterpret_cast< int(*)>(_a[8]))); break;
        case 6: _t->recieveNewShell((*reinterpret_cast< std::string(*)>(_a[1])),(*reinterpret_cast< std::string(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3])),(*reinterpret_cast< float(*)>(_a[4])),(*reinterpret_cast< float(*)>(_a[5])),(*reinterpret_cast< double(*)>(_a[6]))); break;
        case 7: _t->recieveNewBrush((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])),(*reinterpret_cast< bool(*)>(_a[4]))); break;
        case 8: _t->updateBeadTypeTable(); break;
        case 9: _t->updateBeadTable(); break;
        case 10: { double _r = _t->updateM((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< double*>(_a[0]) = std::move(_r); }  break;
        case 11: { double _r = _t->updateTheta((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])),(*reinterpret_cast< double(*)>(_a[3])));
            if (_a[0]) *reinterpret_cast< double*>(_a[0]) = std::move(_r); }  break;
        case 12: _t->on_actionQuit_triggered(); break;
        case 13: _t->on_findUADir_clicked(); break;
        case 14: _t->on_saveNPButton_clicked(); break;
        case 15: _t->on_actionTips_triggered(); break;
        case 16: _t->updateGraphicsWindow(); break;
        case 17: _t->on_updateTables_clicked(); break;
        case 18: _t->on_beadTypeTable_cellChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 19: _t->on_recenterNP_clicked(); break;
        case 20: _t->updateBindingRadii(); break;
        case 21: _t->on_autoBounds_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 22: _t->on_lowerBoundLine_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 23: _t->on_upperBoundLine_textEdited((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 24: _t->on_actionNew_triggered(); break;
        case 25: _t->on_actionLoad_triggered(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject MainWindow::staticMetaObject = { {
    QMetaObject::SuperData::link<QMainWindow::staticMetaObject>(),
    qt_meta_stringdata_MainWindow.data,
    qt_meta_data_MainWindow,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 26)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 26;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 26)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 26;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
