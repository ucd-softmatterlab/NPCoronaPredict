#ifndef CLICKABLESCENE_H
#define CLICKABLESCENE_H

#include <QObject>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>


class ClickableScene : public QGraphicsScene
{
    Q_OBJECT
public:
    explicit ClickableScene(QObject *parent = nullptr);
    virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent * mouseEvent);
signals:
    void sendMouseClickPos( QPointF scenePosLoc);
};

#endif // CLICKABLESCENE_H
