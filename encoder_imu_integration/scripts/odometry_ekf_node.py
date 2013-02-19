#!/usr/bin/env python

"""
This ROS node is responsible taking in wheel encoder
and IMU measurements and outputting an odometry message
with the estimated position combining the two sources.

"""

import roslib; roslib.load_manifest('encoder_imu_integration')
import rospy
from std_msgs.msg import String


def talker():
    pub = rospy.Publisher('chatter', String)
    rospy.init_node('talker')
    while not rospy.is_shutdown():
        str = "hello world %s"%rospy.get_time()
        rospy.loginfo(str)
        pub.publish(String(str))
        rospy.sleep(1.0)
if __name__ == '__main__':
    try:
        talker()
    except rospy.ROSInterruptException: pass



#!/usr/bin/env python
import roslib; roslib.load_manifest('beginner_tutorials')
import rospy
from std_msgs.msg import String
def callback(data):
    rospy.loginfo(rospy.get_name()+"I heard %s",data.data)

def listener():
    rospy.init_node('listener', anonymous=True)
    rospy.Subscriber("chatter", String, callback)
    rospy.spin()

if __name__ == '__main__':
    listener()