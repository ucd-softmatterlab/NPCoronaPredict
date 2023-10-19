#include "clickablescene.h"
#include <QDebug>

ClickableScene::ClickableScene(QObject *parent)
    : QGraphicsScene{parent}
{

}
void ClickableScene::mouseReleaseEvent(QGraphicsSceneMouseEvent * me){
  emit(  sendMouseClickPos(  me->scenePos() ) );
}

