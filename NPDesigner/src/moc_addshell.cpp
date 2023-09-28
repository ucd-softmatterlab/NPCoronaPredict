/****************************************************************************
** Meta object code from reading C++ file 'addshell.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "addshell.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'addshell.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_AddShell_t {
    QByteArrayData data[13];
    char stringdata0[179];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_AddShell_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_AddShell_t qt_meta_stringdata_AddShell = {
    {
QT_MOC_LITERAL(0, 0, 8), // "AddShell"
QT_MOC_LITERAL(1, 9, 12), // "sendNewShell"
QT_MOC_LITERAL(2, 22, 0), // ""
QT_MOC_LITERAL(3, 23, 11), // "std::string"
QT_MOC_LITERAL(4, 35, 13), // "hamakerFileIn"
QT_MOC_LITERAL(5, 49, 12), // "surfaceDirIn"
QT_MOC_LITERAL(6, 62, 11), // "innerRadius"
QT_MOC_LITERAL(7, 74, 11), // "outerRadius"
QT_MOC_LITERAL(8, 86, 18), // "surfacePotentialIn"
QT_MOC_LITERAL(9, 105, 8), // "ljCutoff"
QT_MOC_LITERAL(10, 114, 21), // "on_buttonBox_accepted"
QT_MOC_LITERAL(11, 136, 22), // "on_findHamaker_clicked"
QT_MOC_LITERAL(12, 159, 19) // "on_findSurf_clicked"

    },
    "AddShell\0sendNewShell\0\0std::string\0"
    "hamakerFileIn\0surfaceDirIn\0innerRadius\0"
    "outerRadius\0surfacePotentialIn\0ljCutoff\0"
    "on_buttonBox_accepted\0on_findHamaker_clicked\0"
    "on_findSurf_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_AddShell[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    6,   34,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      10,    0,   47,    2, 0x08 /* Private */,
      11,    0,   48,    2, 0x08 /* Private */,
      12,    0,   49,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Double,    4,    5,    6,    7,    8,    9,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void AddShell::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<AddShell *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->sendNewShell((*reinterpret_cast< std::string(*)>(_a[1])),(*reinterpret_cast< std::string(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3])),(*reinterpret_cast< float(*)>(_a[4])),(*reinterpret_cast< float(*)>(_a[5])),(*reinterpret_cast< double(*)>(_a[6]))); break;
        case 1: _t->on_buttonBox_accepted(); break;
        case 2: _t->on_findHamaker_clicked(); break;
        case 3: _t->on_findSurf_clicked(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (AddShell::*)(std::string , std::string , float , float , float , double );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&AddShell::sendNewShell)) {
                *result = 0;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject AddShell::staticMetaObject = { {
    QMetaObject::SuperData::link<QDialog::staticMetaObject>(),
    qt_meta_stringdata_AddShell.data,
    qt_meta_data_AddShell,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *AddShell::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *AddShell::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_AddShell.stringdata0))
        return static_cast<void*>(this);
    return QDialog::qt_metacast(_clname);
}

int AddShell::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void AddShell::sendNewShell(std::string _t1, std::string _t2, float _t3, float _t4, float _t5, double _t6)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t4))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t5))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t6))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
