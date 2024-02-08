/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../src/mainwindow.h"
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
    QByteArrayData data[65];
    char stringdata0[1257];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 24), // "on_loadUAMButton_clicked"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 17), // "updateHeatmapPlot"
QT_MOC_LITERAL(4, 55, 23), // "on_findUAButton_clicked"
QT_MOC_LITERAL(5, 79, 29), // "on_resultFolderButton_clicked"
QT_MOC_LITERAL(6, 109, 29), // "on_loadMaterialButton_clicked"
QT_MOC_LITERAL(7, 139, 26), // "on_pdbTargetButton_clicked"
QT_MOC_LITERAL(8, 166, 22), // "on_runUAButton_clicked"
QT_MOC_LITERAL(9, 189, 11), // "updateUABox"
QT_MOC_LITERAL(10, 201, 11), // "uaDoneAlert"
QT_MOC_LITERAL(11, 213, 8), // "exitCode"
QT_MOC_LITERAL(12, 222, 20), // "QProcess::ExitStatus"
QT_MOC_LITERAL(13, 243, 10), // "exitStatus"
QT_MOC_LITERAL(14, 254, 15), // "updateEnergyBox"
QT_MOC_LITERAL(15, 270, 27), // "on_phiInputBox_valueChanged"
QT_MOC_LITERAL(16, 298, 4), // "arg1"
QT_MOC_LITERAL(17, 303, 18), // "getSceneMouseClick"
QT_MOC_LITERAL(18, 322, 11), // "scenePosLoc"
QT_MOC_LITERAL(19, 334, 29), // "on_thetaInputBox_valueChanged"
QT_MOC_LITERAL(20, 364, 10), // "phiToIndex"
QT_MOC_LITERAL(21, 375, 6), // "phiVal"
QT_MOC_LITERAL(22, 382, 8), // "deltaVal"
QT_MOC_LITERAL(23, 391, 24), // "on_loadPDBButton_clicked"
QT_MOC_LITERAL(24, 416, 17), // "updateMoleculeBox"
QT_MOC_LITERAL(25, 434, 29), // "on_radiusSpinBox_valueChanged"
QT_MOC_LITERAL(26, 464, 28), // "on_npViewRadius_valueChanged"
QT_MOC_LITERAL(27, 493, 25), // "on_omegaDial_valueChanged"
QT_MOC_LITERAL(28, 519, 5), // "value"
QT_MOC_LITERAL(29, 525, 17), // "checkForMaterials"
QT_MOC_LITERAL(30, 543, 13), // "loadMaterials"
QT_MOC_LITERAL(31, 557, 12), // "materialFile"
QT_MOC_LITERAL(32, 570, 17), // "calcBeadDistances"
QT_MOC_LITERAL(33, 588, 8), // "doRotate"
QT_MOC_LITERAL(34, 597, 22), // "on_colourBoltz_clicked"
QT_MOC_LITERAL(35, 620, 30), // "on_findMinEnergyButton_clicked"
QT_MOC_LITERAL(36, 651, 29), // "on_findBoltzMinButton_clicked"
QT_MOC_LITERAL(37, 681, 23), // "on_colourEnergy_clicked"
QT_MOC_LITERAL(38, 705, 25), // "on_autoNPBox_stateChanged"
QT_MOC_LITERAL(39, 731, 25), // "on_npTargetButton_clicked"
QT_MOC_LITERAL(40, 757, 27), // "on_npcpModeBox_stateChanged"
QT_MOC_LITERAL(41, 785, 45), // "on_mediumEditTable_customCont..."
QT_MOC_LITERAL(42, 831, 3), // "pos"
QT_MOC_LITERAL(43, 835, 19), // "addMoleculeToMedium"
QT_MOC_LITERAL(44, 855, 24), // "removeMoleculeFromMedium"
QT_MOC_LITERAL(45, 880, 26), // "on_mediumNewButton_clicked"
QT_MOC_LITERAL(46, 907, 27), // "on_mediumSaveButton_clicked"
QT_MOC_LITERAL(47, 935, 27), // "on_mediumLoadButton_clicked"
QT_MOC_LITERAL(48, 963, 26), // "on_cancelRunButton_clicked"
QT_MOC_LITERAL(49, 990, 31), // "on_checkStructureButton_clicked"
QT_MOC_LITERAL(50, 1022, 16), // "colourStructures"
QT_MOC_LITERAL(51, 1039, 18), // "colourStructureRow"
QT_MOC_LITERAL(52, 1058, 3), // "row"
QT_MOC_LITERAL(53, 1062, 13), // "startDownload"
QT_MOC_LITERAL(54, 1076, 9), // "targetURL"
QT_MOC_LITERAL(55, 1086, 21), // "downloadReplyFinished"
QT_MOC_LITERAL(56, 1108, 14), // "QNetworkReply*"
QT_MOC_LITERAL(57, 1123, 5), // "reply"
QT_MOC_LITERAL(58, 1129, 30), // "on_mediumEditTable_cellChanged"
QT_MOC_LITERAL(59, 1160, 6), // "column"
QT_MOC_LITERAL(60, 1167, 37), // "on_mediumEditTable_currentCel..."
QT_MOC_LITERAL(61, 1205, 10), // "currentRow"
QT_MOC_LITERAL(62, 1216, 13), // "currentColumn"
QT_MOC_LITERAL(63, 1230, 11), // "previousRow"
QT_MOC_LITERAL(64, 1242, 14) // "previousColumn"

    },
    "MainWindow\0on_loadUAMButton_clicked\0"
    "\0updateHeatmapPlot\0on_findUAButton_clicked\0"
    "on_resultFolderButton_clicked\0"
    "on_loadMaterialButton_clicked\0"
    "on_pdbTargetButton_clicked\0"
    "on_runUAButton_clicked\0updateUABox\0"
    "uaDoneAlert\0exitCode\0QProcess::ExitStatus\0"
    "exitStatus\0updateEnergyBox\0"
    "on_phiInputBox_valueChanged\0arg1\0"
    "getSceneMouseClick\0scenePosLoc\0"
    "on_thetaInputBox_valueChanged\0phiToIndex\0"
    "phiVal\0deltaVal\0on_loadPDBButton_clicked\0"
    "updateMoleculeBox\0on_radiusSpinBox_valueChanged\0"
    "on_npViewRadius_valueChanged\0"
    "on_omegaDial_valueChanged\0value\0"
    "checkForMaterials\0loadMaterials\0"
    "materialFile\0calcBeadDistances\0doRotate\0"
    "on_colourBoltz_clicked\0"
    "on_findMinEnergyButton_clicked\0"
    "on_findBoltzMinButton_clicked\0"
    "on_colourEnergy_clicked\0"
    "on_autoNPBox_stateChanged\0"
    "on_npTargetButton_clicked\0"
    "on_npcpModeBox_stateChanged\0"
    "on_mediumEditTable_customContextMenuRequested\0"
    "pos\0addMoleculeToMedium\0"
    "removeMoleculeFromMedium\0"
    "on_mediumNewButton_clicked\0"
    "on_mediumSaveButton_clicked\0"
    "on_mediumLoadButton_clicked\0"
    "on_cancelRunButton_clicked\0"
    "on_checkStructureButton_clicked\0"
    "colourStructures\0colourStructureRow\0"
    "row\0startDownload\0targetURL\0"
    "downloadReplyFinished\0QNetworkReply*\0"
    "reply\0on_mediumEditTable_cellChanged\0"
    "column\0on_mediumEditTable_currentCellChanged\0"
    "currentRow\0currentColumn\0previousRow\0"
    "previousColumn"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      43,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  229,    2, 0x08 /* Private */,
       3,    0,  230,    2, 0x08 /* Private */,
       4,    0,  231,    2, 0x08 /* Private */,
       5,    0,  232,    2, 0x08 /* Private */,
       6,    0,  233,    2, 0x08 /* Private */,
       7,    0,  234,    2, 0x08 /* Private */,
       8,    0,  235,    2, 0x08 /* Private */,
       9,    0,  236,    2, 0x08 /* Private */,
      10,    2,  237,    2, 0x08 /* Private */,
      14,    0,  242,    2, 0x08 /* Private */,
      15,    1,  243,    2, 0x08 /* Private */,
      17,    1,  246,    2, 0x08 /* Private */,
      19,    1,  249,    2, 0x08 /* Private */,
      20,    2,  252,    2, 0x08 /* Private */,
      23,    0,  257,    2, 0x08 /* Private */,
      24,    0,  258,    2, 0x08 /* Private */,
      25,    1,  259,    2, 0x08 /* Private */,
      26,    1,  262,    2, 0x08 /* Private */,
      27,    1,  265,    2, 0x08 /* Private */,
      29,    0,  268,    2, 0x08 /* Private */,
      30,    1,  269,    2, 0x08 /* Private */,
      32,    1,  272,    2, 0x08 /* Private */,
      34,    0,  275,    2, 0x08 /* Private */,
      35,    0,  276,    2, 0x08 /* Private */,
      36,    0,  277,    2, 0x08 /* Private */,
      37,    0,  278,    2, 0x08 /* Private */,
      38,    1,  279,    2, 0x08 /* Private */,
      39,    0,  282,    2, 0x08 /* Private */,
      40,    1,  283,    2, 0x08 /* Private */,
      41,    1,  286,    2, 0x08 /* Private */,
      43,    0,  289,    2, 0x08 /* Private */,
      44,    0,  290,    2, 0x08 /* Private */,
      45,    0,  291,    2, 0x08 /* Private */,
      46,    0,  292,    2, 0x08 /* Private */,
      47,    0,  293,    2, 0x08 /* Private */,
      48,    0,  294,    2, 0x08 /* Private */,
      49,    0,  295,    2, 0x08 /* Private */,
      50,    0,  296,    2, 0x08 /* Private */,
      51,    1,  297,    2, 0x08 /* Private */,
      53,    1,  300,    2, 0x08 /* Private */,
      55,    1,  303,    2, 0x08 /* Private */,
      58,    2,  306,    2, 0x08 /* Private */,
      60,    4,  311,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, 0x80000000 | 12,   11,   13,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::QPointF,   18,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Int, QMetaType::Double, QMetaType::Double,   21,   22,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::Int,   28,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,   31,
    QMetaType::Void, QMetaType::Bool,   33,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::QPoint,   42,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   52,
    QMetaType::Void, QMetaType::QUrl,   54,
    QMetaType::Void, 0x80000000 | 56,   57,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   52,   59,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::Int, QMetaType::Int,   61,   62,   63,   64,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MainWindow *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->on_loadUAMButton_clicked(); break;
        case 1: _t->updateHeatmapPlot(); break;
        case 2: _t->on_findUAButton_clicked(); break;
        case 3: _t->on_resultFolderButton_clicked(); break;
        case 4: _t->on_loadMaterialButton_clicked(); break;
        case 5: _t->on_pdbTargetButton_clicked(); break;
        case 6: _t->on_runUAButton_clicked(); break;
        case 7: _t->updateUABox(); break;
        case 8: _t->uaDoneAlert((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QProcess::ExitStatus(*)>(_a[2]))); break;
        case 9: _t->updateEnergyBox(); break;
        case 10: _t->on_phiInputBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->getSceneMouseClick((*reinterpret_cast< QPointF(*)>(_a[1]))); break;
        case 12: _t->on_thetaInputBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: { int _r = _t->phiToIndex((*reinterpret_cast< double(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 14: _t->on_loadPDBButton_clicked(); break;
        case 15: _t->updateMoleculeBox(); break;
        case 16: _t->on_radiusSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: _t->on_npViewRadius_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 18: _t->on_omegaDial_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 19: _t->checkForMaterials(); break;
        case 20: _t->loadMaterials((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 21: _t->calcBeadDistances((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 22: _t->on_colourBoltz_clicked(); break;
        case 23: _t->on_findMinEnergyButton_clicked(); break;
        case 24: _t->on_findBoltzMinButton_clicked(); break;
        case 25: _t->on_colourEnergy_clicked(); break;
        case 26: _t->on_autoNPBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 27: _t->on_npTargetButton_clicked(); break;
        case 28: _t->on_npcpModeBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 29: _t->on_mediumEditTable_customContextMenuRequested((*reinterpret_cast< const QPoint(*)>(_a[1]))); break;
        case 30: _t->addMoleculeToMedium(); break;
        case 31: _t->removeMoleculeFromMedium(); break;
        case 32: _t->on_mediumNewButton_clicked(); break;
        case 33: _t->on_mediumSaveButton_clicked(); break;
        case 34: _t->on_mediumLoadButton_clicked(); break;
        case 35: _t->on_cancelRunButton_clicked(); break;
        case 36: _t->on_checkStructureButton_clicked(); break;
        case 37: _t->colourStructures(); break;
        case 38: _t->colourStructureRow((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 39: _t->startDownload((*reinterpret_cast< QUrl(*)>(_a[1]))); break;
        case 40: _t->downloadReplyFinished((*reinterpret_cast< QNetworkReply*(*)>(_a[1]))); break;
        case 41: _t->on_mediumEditTable_cellChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 42: _t->on_mediumEditTable_currentCellChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 40:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QNetworkReply* >(); break;
            }
            break;
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
        if (_id < 43)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 43;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 43)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 43;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
