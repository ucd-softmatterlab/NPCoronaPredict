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
    QByteArrayData data[79];
    char stringdata0[1594];
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
QT_MOC_LITERAL(32, 570, 7), // "doAlert"
QT_MOC_LITERAL(33, 578, 17), // "calcBeadDistances"
QT_MOC_LITERAL(34, 596, 8), // "doRotate"
QT_MOC_LITERAL(35, 605, 22), // "calcBeadBoltzDistances"
QT_MOC_LITERAL(36, 628, 22), // "on_colourBoltz_clicked"
QT_MOC_LITERAL(37, 651, 30), // "on_findMinEnergyButton_clicked"
QT_MOC_LITERAL(38, 682, 29), // "on_findBoltzMinButton_clicked"
QT_MOC_LITERAL(39, 712, 23), // "on_colourEnergy_clicked"
QT_MOC_LITERAL(40, 736, 25), // "on_autoNPBox_stateChanged"
QT_MOC_LITERAL(41, 762, 25), // "on_npTargetButton_clicked"
QT_MOC_LITERAL(42, 788, 27), // "on_npcpModeBox_stateChanged"
QT_MOC_LITERAL(43, 816, 45), // "on_mediumEditTable_customCont..."
QT_MOC_LITERAL(44, 862, 3), // "pos"
QT_MOC_LITERAL(45, 866, 19), // "addMoleculeToMedium"
QT_MOC_LITERAL(46, 886, 24), // "removeMoleculeFromMedium"
QT_MOC_LITERAL(47, 911, 26), // "on_mediumNewButton_clicked"
QT_MOC_LITERAL(48, 938, 27), // "on_mediumSaveButton_clicked"
QT_MOC_LITERAL(49, 966, 27), // "on_mediumLoadButton_clicked"
QT_MOC_LITERAL(50, 994, 26), // "on_cancelRunButton_clicked"
QT_MOC_LITERAL(51, 1021, 31), // "on_checkStructureButton_clicked"
QT_MOC_LITERAL(52, 1053, 16), // "colourStructures"
QT_MOC_LITERAL(53, 1070, 18), // "colourStructureRow"
QT_MOC_LITERAL(54, 1089, 3), // "row"
QT_MOC_LITERAL(55, 1093, 13), // "startDownload"
QT_MOC_LITERAL(56, 1107, 9), // "targetURL"
QT_MOC_LITERAL(57, 1117, 21), // "downloadReplyFinished"
QT_MOC_LITERAL(58, 1139, 14), // "QNetworkReply*"
QT_MOC_LITERAL(59, 1154, 5), // "reply"
QT_MOC_LITERAL(60, 1160, 30), // "on_mediumEditTable_cellChanged"
QT_MOC_LITERAL(61, 1191, 6), // "column"
QT_MOC_LITERAL(62, 1198, 37), // "on_mediumEditTable_currentCel..."
QT_MOC_LITERAL(63, 1236, 10), // "currentRow"
QT_MOC_LITERAL(64, 1247, 13), // "currentColumn"
QT_MOC_LITERAL(65, 1261, 11), // "previousRow"
QT_MOC_LITERAL(66, 1273, 14), // "previousColumn"
QT_MOC_LITERAL(67, 1288, 26), // "on_showCOMBox_stateChanged"
QT_MOC_LITERAL(68, 1315, 28), // "on_opacitySlider_sliderMoved"
QT_MOC_LITERAL(69, 1344, 8), // "position"
QT_MOC_LITERAL(70, 1353, 29), // "on_opacitySlider_valueChanged"
QT_MOC_LITERAL(71, 1383, 28), // "on_loadBeadmapButton_clicked"
QT_MOC_LITERAL(72, 1412, 15), // "loadBeadSetFile"
QT_MOC_LITERAL(73, 1428, 10), // "targetFile"
QT_MOC_LITERAL(74, 1439, 29), // "on_showChargeBox_stateChanged"
QT_MOC_LITERAL(75, 1469, 33), // "on_npShapeBox_currentIndexCha..."
QT_MOC_LITERAL(76, 1503, 5), // "index"
QT_MOC_LITERAL(77, 1509, 39), // "on_materialDropdown_currentIn..."
QT_MOC_LITERAL(78, 1549, 44) // "on_npTargetShapeOverride_curr..."

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
    "materialFile\0doAlert\0calcBeadDistances\0"
    "doRotate\0calcBeadBoltzDistances\0"
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
    "previousColumn\0on_showCOMBox_stateChanged\0"
    "on_opacitySlider_sliderMoved\0position\0"
    "on_opacitySlider_valueChanged\0"
    "on_loadBeadmapButton_clicked\0"
    "loadBeadSetFile\0targetFile\0"
    "on_showChargeBox_stateChanged\0"
    "on_npShapeBox_currentIndexChanged\0"
    "index\0on_materialDropdown_currentIndexChanged\0"
    "on_npTargetShapeOverride_currentIndexChanged"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      53,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  279,    2, 0x08 /* Private */,
       3,    0,  280,    2, 0x08 /* Private */,
       4,    0,  281,    2, 0x08 /* Private */,
       5,    0,  282,    2, 0x08 /* Private */,
       6,    0,  283,    2, 0x08 /* Private */,
       7,    0,  284,    2, 0x08 /* Private */,
       8,    0,  285,    2, 0x08 /* Private */,
       9,    0,  286,    2, 0x08 /* Private */,
      10,    2,  287,    2, 0x08 /* Private */,
      14,    0,  292,    2, 0x08 /* Private */,
      15,    1,  293,    2, 0x08 /* Private */,
      17,    1,  296,    2, 0x08 /* Private */,
      19,    1,  299,    2, 0x08 /* Private */,
      20,    2,  302,    2, 0x08 /* Private */,
      23,    0,  307,    2, 0x08 /* Private */,
      24,    0,  308,    2, 0x08 /* Private */,
      25,    1,  309,    2, 0x08 /* Private */,
      26,    1,  312,    2, 0x08 /* Private */,
      27,    1,  315,    2, 0x08 /* Private */,
      29,    0,  318,    2, 0x08 /* Private */,
      30,    2,  319,    2, 0x08 /* Private */,
      33,    1,  324,    2, 0x08 /* Private */,
      35,    1,  327,    2, 0x08 /* Private */,
      36,    0,  330,    2, 0x08 /* Private */,
      37,    0,  331,    2, 0x08 /* Private */,
      38,    0,  332,    2, 0x08 /* Private */,
      39,    0,  333,    2, 0x08 /* Private */,
      40,    1,  334,    2, 0x08 /* Private */,
      41,    0,  337,    2, 0x08 /* Private */,
      42,    1,  338,    2, 0x08 /* Private */,
      43,    1,  341,    2, 0x08 /* Private */,
      45,    0,  344,    2, 0x08 /* Private */,
      46,    0,  345,    2, 0x08 /* Private */,
      47,    0,  346,    2, 0x08 /* Private */,
      48,    0,  347,    2, 0x08 /* Private */,
      49,    0,  348,    2, 0x08 /* Private */,
      50,    0,  349,    2, 0x08 /* Private */,
      51,    0,  350,    2, 0x08 /* Private */,
      52,    0,  351,    2, 0x08 /* Private */,
      53,    1,  352,    2, 0x08 /* Private */,
      55,    1,  355,    2, 0x08 /* Private */,
      57,    1,  358,    2, 0x08 /* Private */,
      60,    2,  361,    2, 0x08 /* Private */,
      62,    4,  366,    2, 0x08 /* Private */,
      67,    1,  375,    2, 0x08 /* Private */,
      68,    1,  378,    2, 0x08 /* Private */,
      70,    1,  381,    2, 0x08 /* Private */,
      71,    0,  384,    2, 0x08 /* Private */,
      72,    1,  385,    2, 0x08 /* Private */,
      74,    1,  388,    2, 0x08 /* Private */,
      75,    1,  391,    2, 0x08 /* Private */,
      77,    1,  394,    2, 0x08 /* Private */,
      78,    1,  397,    2, 0x08 /* Private */,

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
    QMetaType::Void, QMetaType::QString, QMetaType::Bool,   31,   32,
    QMetaType::Void, QMetaType::Bool,   34,
    QMetaType::Void, QMetaType::Bool,   34,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::QPoint,   44,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   54,
    QMetaType::Void, QMetaType::QUrl,   56,
    QMetaType::Void, 0x80000000 | 58,   59,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   54,   61,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::Int, QMetaType::Int,   63,   64,   65,   66,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::Int,   69,
    QMetaType::Void, QMetaType::Int,   28,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,   73,
    QMetaType::Void, QMetaType::Int,   16,
    QMetaType::Void, QMetaType::Int,   76,
    QMetaType::Void, QMetaType::Int,   76,
    QMetaType::Void, QMetaType::Int,   76,

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
        case 20: _t->loadMaterials((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 21: _t->calcBeadDistances((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 22: _t->calcBeadBoltzDistances((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 23: _t->on_colourBoltz_clicked(); break;
        case 24: _t->on_findMinEnergyButton_clicked(); break;
        case 25: _t->on_findBoltzMinButton_clicked(); break;
        case 26: _t->on_colourEnergy_clicked(); break;
        case 27: _t->on_autoNPBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 28: _t->on_npTargetButton_clicked(); break;
        case 29: _t->on_npcpModeBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 30: _t->on_mediumEditTable_customContextMenuRequested((*reinterpret_cast< const QPoint(*)>(_a[1]))); break;
        case 31: _t->addMoleculeToMedium(); break;
        case 32: _t->removeMoleculeFromMedium(); break;
        case 33: _t->on_mediumNewButton_clicked(); break;
        case 34: _t->on_mediumSaveButton_clicked(); break;
        case 35: _t->on_mediumLoadButton_clicked(); break;
        case 36: _t->on_cancelRunButton_clicked(); break;
        case 37: _t->on_checkStructureButton_clicked(); break;
        case 38: _t->colourStructures(); break;
        case 39: _t->colourStructureRow((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 40: _t->startDownload((*reinterpret_cast< QUrl(*)>(_a[1]))); break;
        case 41: _t->downloadReplyFinished((*reinterpret_cast< QNetworkReply*(*)>(_a[1]))); break;
        case 42: _t->on_mediumEditTable_cellChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 43: _t->on_mediumEditTable_currentCellChanged((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4]))); break;
        case 44: _t->on_showCOMBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 45: _t->on_opacitySlider_sliderMoved((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 46: _t->on_opacitySlider_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 47: _t->on_loadBeadmapButton_clicked(); break;
        case 48: _t->loadBeadSetFile((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 49: _t->on_showChargeBox_stateChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 50: _t->on_npShapeBox_currentIndexChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 51: _t->on_materialDropdown_currentIndexChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 52: _t->on_npTargetShapeOverride_currentIndexChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 41:
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
        if (_id < 53)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 53;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 53)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 53;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
