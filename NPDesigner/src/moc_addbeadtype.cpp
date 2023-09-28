/****************************************************************************
** Meta object code from reading C++ file 'addbeadtype.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.15.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "addbeadtype.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'addbeadtype.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.15.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_AddBeadType_t {
    QByteArrayData data[15];
    char stringdata0[232];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_AddBeadType_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_AddBeadType_t qt_meta_stringdata_AddBeadType = {
    {
QT_MOC_LITERAL(0, 0, 11), // "AddBeadType"
QT_MOC_LITERAL(1, 12, 17), // "sendNewNPBeadType"
QT_MOC_LITERAL(2, 30, 0), // ""
QT_MOC_LITERAL(3, 31, 11), // "std::string"
QT_MOC_LITERAL(4, 43, 13), // "hamakerFileIn"
QT_MOC_LITERAL(5, 57, 12), // "surfaceDirIn"
QT_MOC_LITERAL(6, 70, 8), // "radiusIn"
QT_MOC_LITERAL(7, 79, 18), // "surfacePotentialIn"
QT_MOC_LITERAL(8, 98, 12), // "surfFactorIn"
QT_MOC_LITERAL(9, 111, 12), // "coreFactorIn"
QT_MOC_LITERAL(10, 124, 10), // "ljCutoffIn"
QT_MOC_LITERAL(11, 135, 20), // "correctionOverrideIn"
QT_MOC_LITERAL(12, 156, 21), // "on_buttonBox_accepted"
QT_MOC_LITERAL(13, 178, 28), // "on_surfFindDirButton_clicked"
QT_MOC_LITERAL(14, 207, 24) // "on_hamFindButton_clicked"

    },
    "AddBeadType\0sendNewNPBeadType\0\0"
    "std::string\0hamakerFileIn\0surfaceDirIn\0"
    "radiusIn\0surfacePotentialIn\0surfFactorIn\0"
    "coreFactorIn\0ljCutoffIn\0correctionOverrideIn\0"
    "on_buttonBox_accepted\0"
    "on_surfFindDirButton_clicked\0"
    "on_hamFindButton_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_AddBeadType[] = {

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
       1,    8,   34,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      12,    0,   51,    2, 0x08 /* Private */,
      13,    0,   52,    2, 0x08 /* Private */,
      14,    0,   53,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Float, QMetaType::Int,    4,    5,    6,    7,    8,    9,   10,   11,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void AddBeadType::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<AddBeadType *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->sendNewNPBeadType((*reinterpret_cast< std::string(*)>(_a[1])),(*reinterpret_cast< std::string(*)>(_a[2])),(*reinterpret_cast< float(*)>(_a[3])),(*reinterpret_cast< float(*)>(_a[4])),(*reinterpret_cast< float(*)>(_a[5])),(*reinterpret_cast< float(*)>(_a[6])),(*reinterpret_cast< float(*)>(_a[7])),(*reinterpret_cast< int(*)>(_a[8]))); break;
        case 1: _t->on_buttonBox_accepted(); break;
        case 2: _t->on_surfFindDirButton_clicked(); break;
        case 3: _t->on_hamFindButton_clicked(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (AddBeadType::*)(std::string , std::string , float , float , float , float , float , int );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&AddBeadType::sendNewNPBeadType)) {
                *result = 0;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject AddBeadType::staticMetaObject = { {
    QMetaObject::SuperData::link<QDialog::staticMetaObject>(),
    qt_meta_stringdata_AddBeadType.data,
    qt_meta_data_AddBeadType,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *AddBeadType::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *AddBeadType::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_AddBeadType.stringdata0))
        return static_cast<void*>(this);
    return QDialog::qt_metacast(_clname);
}

int AddBeadType::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
void AddBeadType::sendNewNPBeadType(std::string _t1, std::string _t2, float _t3, float _t4, float _t5, float _t6, float _t7, int _t8)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t3))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t4))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t5))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t6))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t7))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t8))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
