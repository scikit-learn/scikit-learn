import sys
from _typeshed import ReadableBuffer, Self, WriteableBuffer
from collections.abc import Iterable
from enum import IntEnum, IntFlag
from io import RawIOBase
from typing import Any, BinaryIO, TextIO, TypeVar, overload
from typing_extensions import Literal

# Ideally, we'd just do "from _socket import *". Unfortunately, socket
# overrides some definitions from _socket incompatibly. mypy incorrectly
# prefers the definitions from _socket over those defined here.
import _socket
from _socket import (
    _FD,
    CMSG_LEN as CMSG_LEN,
    CMSG_SPACE as CMSG_SPACE,
    EAI_ADDRFAMILY as EAI_ADDRFAMILY,
    EAI_AGAIN as EAI_AGAIN,
    EAI_BADFLAGS as EAI_BADFLAGS,
    EAI_BADHINTS as EAI_BADHINTS,
    EAI_FAIL as EAI_FAIL,
    EAI_FAMILY as EAI_FAMILY,
    EAI_MAX as EAI_MAX,
    EAI_MEMORY as EAI_MEMORY,
    EAI_NODATA as EAI_NODATA,
    EAI_NONAME as EAI_NONAME,
    EAI_OVERFLOW as EAI_OVERFLOW,
    EAI_PROTOCOL as EAI_PROTOCOL,
    EAI_SERVICE as EAI_SERVICE,
    EAI_SOCKTYPE as EAI_SOCKTYPE,
    EAI_SYSTEM as EAI_SYSTEM,
    INADDR_ALLHOSTS_GROUP as INADDR_ALLHOSTS_GROUP,
    INADDR_ANY as INADDR_ANY,
    INADDR_BROADCAST as INADDR_BROADCAST,
    INADDR_LOOPBACK as INADDR_LOOPBACK,
    INADDR_MAX_LOCAL_GROUP as INADDR_MAX_LOCAL_GROUP,
    INADDR_NONE as INADDR_NONE,
    INADDR_UNSPEC_GROUP as INADDR_UNSPEC_GROUP,
    IP_ADD_MEMBERSHIP as IP_ADD_MEMBERSHIP,
    IP_DEFAULT_MULTICAST_LOOP as IP_DEFAULT_MULTICAST_LOOP,
    IP_DEFAULT_MULTICAST_TTL as IP_DEFAULT_MULTICAST_TTL,
    IP_DROP_MEMBERSHIP as IP_DROP_MEMBERSHIP,
    IP_HDRINCL as IP_HDRINCL,
    IP_MAX_MEMBERSHIPS as IP_MAX_MEMBERSHIPS,
    IP_MULTICAST_IF as IP_MULTICAST_IF,
    IP_MULTICAST_LOOP as IP_MULTICAST_LOOP,
    IP_MULTICAST_TTL as IP_MULTICAST_TTL,
    IP_OPTIONS as IP_OPTIONS,
    IP_RECVDSTADDR as IP_RECVDSTADDR,
    IP_RECVOPTS as IP_RECVOPTS,
    IP_RECVRETOPTS as IP_RECVRETOPTS,
    IP_RETOPTS as IP_RETOPTS,
    IP_TOS as IP_TOS,
    IP_TRANSPARENT as IP_TRANSPARENT,
    IP_TTL as IP_TTL,
    IPPORT_RESERVED as IPPORT_RESERVED,
    IPPORT_USERRESERVED as IPPORT_USERRESERVED,
    IPPROTO_AH as IPPROTO_AH,
    IPPROTO_BIP as IPPROTO_BIP,
    IPPROTO_DSTOPTS as IPPROTO_DSTOPTS,
    IPPROTO_EGP as IPPROTO_EGP,
    IPPROTO_EON as IPPROTO_EON,
    IPPROTO_ESP as IPPROTO_ESP,
    IPPROTO_FRAGMENT as IPPROTO_FRAGMENT,
    IPPROTO_GGP as IPPROTO_GGP,
    IPPROTO_GRE as IPPROTO_GRE,
    IPPROTO_HELLO as IPPROTO_HELLO,
    IPPROTO_HOPOPTS as IPPROTO_HOPOPTS,
    IPPROTO_ICMP as IPPROTO_ICMP,
    IPPROTO_ICMPV6 as IPPROTO_ICMPV6,
    IPPROTO_IDP as IPPROTO_IDP,
    IPPROTO_IGMP as IPPROTO_IGMP,
    IPPROTO_IP as IPPROTO_IP,
    IPPROTO_IPCOMP as IPPROTO_IPCOMP,
    IPPROTO_IPIP as IPPROTO_IPIP,
    IPPROTO_IPV4 as IPPROTO_IPV4,
    IPPROTO_IPV6 as IPPROTO_IPV6,
    IPPROTO_MAX as IPPROTO_MAX,
    IPPROTO_MOBILE as IPPROTO_MOBILE,
    IPPROTO_ND as IPPROTO_ND,
    IPPROTO_NONE as IPPROTO_NONE,
    IPPROTO_PIM as IPPROTO_PIM,
    IPPROTO_PUP as IPPROTO_PUP,
    IPPROTO_RAW as IPPROTO_RAW,
    IPPROTO_ROUTING as IPPROTO_ROUTING,
    IPPROTO_RSVP as IPPROTO_RSVP,
    IPPROTO_SCTP as IPPROTO_SCTP,
    IPPROTO_TCP as IPPROTO_TCP,
    IPPROTO_TP as IPPROTO_TP,
    IPPROTO_UDP as IPPROTO_UDP,
    IPPROTO_VRRP as IPPROTO_VRRP,
    IPPROTO_XTP as IPPROTO_XTP,
    IPV6_CHECKSUM as IPV6_CHECKSUM,
    IPV6_DONTFRAG as IPV6_DONTFRAG,
    IPV6_DSTOPTS as IPV6_DSTOPTS,
    IPV6_HOPLIMIT as IPV6_HOPLIMIT,
    IPV6_HOPOPTS as IPV6_HOPOPTS,
    IPV6_JOIN_GROUP as IPV6_JOIN_GROUP,
    IPV6_LEAVE_GROUP as IPV6_LEAVE_GROUP,
    IPV6_MULTICAST_HOPS as IPV6_MULTICAST_HOPS,
    IPV6_MULTICAST_IF as IPV6_MULTICAST_IF,
    IPV6_MULTICAST_LOOP as IPV6_MULTICAST_LOOP,
    IPV6_NEXTHOP as IPV6_NEXTHOP,
    IPV6_PATHMTU as IPV6_PATHMTU,
    IPV6_PKTINFO as IPV6_PKTINFO,
    IPV6_RECVDSTOPTS as IPV6_RECVDSTOPTS,
    IPV6_RECVHOPLIMIT as IPV6_RECVHOPLIMIT,
    IPV6_RECVHOPOPTS as IPV6_RECVHOPOPTS,
    IPV6_RECVPATHMTU as IPV6_RECVPATHMTU,
    IPV6_RECVPKTINFO as IPV6_RECVPKTINFO,
    IPV6_RECVRTHDR as IPV6_RECVRTHDR,
    IPV6_RECVTCLASS as IPV6_RECVTCLASS,
    IPV6_RTHDR as IPV6_RTHDR,
    IPV6_RTHDR_TYPE_0 as IPV6_RTHDR_TYPE_0,
    IPV6_RTHDRDSTOPTS as IPV6_RTHDRDSTOPTS,
    IPV6_TCLASS as IPV6_TCLASS,
    IPV6_UNICAST_HOPS as IPV6_UNICAST_HOPS,
    IPV6_USE_MIN_MTU as IPV6_USE_MIN_MTU,
    IPV6_V6ONLY as IPV6_V6ONLY,
    IPX_TYPE as IPX_TYPE,
    LOCAL_PEERCRED as LOCAL_PEERCRED,
    NI_DGRAM as NI_DGRAM,
    NI_MAXHOST as NI_MAXHOST,
    NI_MAXSERV as NI_MAXSERV,
    NI_NAMEREQD as NI_NAMEREQD,
    NI_NOFQDN as NI_NOFQDN,
    NI_NUMERICHOST as NI_NUMERICHOST,
    NI_NUMERICSERV as NI_NUMERICSERV,
    SCM_CREDENTIALS as SCM_CREDENTIALS,
    SCM_CREDS as SCM_CREDS,
    SCM_RIGHTS as SCM_RIGHTS,
    SHUT_RD as SHUT_RD,
    SHUT_RDWR as SHUT_RDWR,
    SHUT_WR as SHUT_WR,
    SO_ACCEPTCONN as SO_ACCEPTCONN,
    SO_BINDTODEVICE as SO_BINDTODEVICE,
    SO_BROADCAST as SO_BROADCAST,
    SO_DEBUG as SO_DEBUG,
    SO_DONTROUTE as SO_DONTROUTE,
    SO_ERROR as SO_ERROR,
    SO_EXCLUSIVEADDRUSE as SO_EXCLUSIVEADDRUSE,
    SO_KEEPALIVE as SO_KEEPALIVE,
    SO_LINGER as SO_LINGER,
    SO_MARK as SO_MARK,
    SO_OOBINLINE as SO_OOBINLINE,
    SO_PASSCRED as SO_PASSCRED,
    SO_PEERCRED as SO_PEERCRED,
    SO_PRIORITY as SO_PRIORITY,
    SO_RCVBUF as SO_RCVBUF,
    SO_RCVLOWAT as SO_RCVLOWAT,
    SO_RCVTIMEO as SO_RCVTIMEO,
    SO_REUSEADDR as SO_REUSEADDR,
    SO_REUSEPORT as SO_REUSEPORT,
    SO_SETFIB as SO_SETFIB,
    SO_SNDBUF as SO_SNDBUF,
    SO_SNDLOWAT as SO_SNDLOWAT,
    SO_SNDTIMEO as SO_SNDTIMEO,
    SO_TYPE as SO_TYPE,
    SO_USELOOPBACK as SO_USELOOPBACK,
    SOL_ATALK as SOL_ATALK,
    SOL_AX25 as SOL_AX25,
    SOL_HCI as SOL_HCI,
    SOL_IP as SOL_IP,
    SOL_IPX as SOL_IPX,
    SOL_NETROM as SOL_NETROM,
    SOL_ROSE as SOL_ROSE,
    SOL_SOCKET as SOL_SOCKET,
    SOL_TCP as SOL_TCP,
    SOL_UDP as SOL_UDP,
    SOMAXCONN as SOMAXCONN,
    TCP_CORK as TCP_CORK,
    TCP_DEFER_ACCEPT as TCP_DEFER_ACCEPT,
    TCP_FASTOPEN as TCP_FASTOPEN,
    TCP_INFO as TCP_INFO,
    TCP_KEEPCNT as TCP_KEEPCNT,
    TCP_KEEPIDLE as TCP_KEEPIDLE,
    TCP_KEEPINTVL as TCP_KEEPINTVL,
    TCP_LINGER2 as TCP_LINGER2,
    TCP_MAXSEG as TCP_MAXSEG,
    TCP_NODELAY as TCP_NODELAY,
    TCP_QUICKACK as TCP_QUICKACK,
    TCP_SYNCNT as TCP_SYNCNT,
    TCP_WINDOW_CLAMP as TCP_WINDOW_CLAMP,
    SocketType as SocketType,
    _Address as _Address,
    _RetAddress as _RetAddress,
    dup as dup,
    error as error,
    gaierror as gaierror,
    getdefaulttimeout as getdefaulttimeout,
    gethostbyaddr as gethostbyaddr,
    gethostbyname as gethostbyname,
    gethostbyname_ex as gethostbyname_ex,
    gethostname as gethostname,
    getnameinfo as getnameinfo,
    getprotobyname as getprotobyname,
    getservbyname as getservbyname,
    getservbyport as getservbyport,
    has_ipv6 as has_ipv6,
    herror as herror,
    htonl as htonl,
    htons as htons,
    inet_aton as inet_aton,
    inet_ntoa as inet_ntoa,
    inet_ntop as inet_ntop,
    inet_pton as inet_pton,
    ntohl as ntohl,
    ntohs as ntohs,
    setdefaulttimeout as setdefaulttimeout,
    timeout as timeout,
)

if sys.version_info >= (3, 7):
    from _socket import close as close
if sys.platform != "win32":
    from _socket import sethostname as sethostname
if sys.platform != "win32" or sys.version_info >= (3, 8):
    from _socket import if_indextoname as if_indextoname, if_nameindex as if_nameindex, if_nametoindex as if_nametoindex
if sys.platform == "linux":
    from _socket import (
        ALG_OP_DECRYPT as ALG_OP_DECRYPT,
        ALG_OP_ENCRYPT as ALG_OP_ENCRYPT,
        ALG_OP_SIGN as ALG_OP_SIGN,
        ALG_OP_VERIFY as ALG_OP_VERIFY,
        ALG_SET_AEAD_ASSOCLEN as ALG_SET_AEAD_ASSOCLEN,
        ALG_SET_AEAD_AUTHSIZE as ALG_SET_AEAD_AUTHSIZE,
        ALG_SET_IV as ALG_SET_IV,
        ALG_SET_KEY as ALG_SET_KEY,
        ALG_SET_OP as ALG_SET_OP,
        ALG_SET_PUBKEY as ALG_SET_PUBKEY,
        CAN_BCM as CAN_BCM,
        CAN_BCM_RX_CHANGED as CAN_BCM_RX_CHANGED,
        CAN_BCM_RX_DELETE as CAN_BCM_RX_DELETE,
        CAN_BCM_RX_READ as CAN_BCM_RX_READ,
        CAN_BCM_RX_SETUP as CAN_BCM_RX_SETUP,
        CAN_BCM_RX_STATUS as CAN_BCM_RX_STATUS,
        CAN_BCM_RX_TIMEOUT as CAN_BCM_RX_TIMEOUT,
        CAN_BCM_TX_DELETE as CAN_BCM_TX_DELETE,
        CAN_BCM_TX_EXPIRED as CAN_BCM_TX_EXPIRED,
        CAN_BCM_TX_READ as CAN_BCM_TX_READ,
        CAN_BCM_TX_SEND as CAN_BCM_TX_SEND,
        CAN_BCM_TX_SETUP as CAN_BCM_TX_SETUP,
        CAN_BCM_TX_STATUS as CAN_BCM_TX_STATUS,
        CAN_EFF_FLAG as CAN_EFF_FLAG,
        CAN_EFF_MASK as CAN_EFF_MASK,
        CAN_ERR_FLAG as CAN_ERR_FLAG,
        CAN_ERR_MASK as CAN_ERR_MASK,
        CAN_RAW as CAN_RAW,
        CAN_RAW_ERR_FILTER as CAN_RAW_ERR_FILTER,
        CAN_RAW_FD_FRAMES as CAN_RAW_FD_FRAMES,
        CAN_RAW_FILTER as CAN_RAW_FILTER,
        CAN_RAW_LOOPBACK as CAN_RAW_LOOPBACK,
        CAN_RAW_RECV_OWN_MSGS as CAN_RAW_RECV_OWN_MSGS,
        CAN_RTR_FLAG as CAN_RTR_FLAG,
        CAN_SFF_MASK as CAN_SFF_MASK,
        PACKET_BROADCAST as PACKET_BROADCAST,
        PACKET_FASTROUTE as PACKET_FASTROUTE,
        PACKET_HOST as PACKET_HOST,
        PACKET_LOOPBACK as PACKET_LOOPBACK,
        PACKET_MULTICAST as PACKET_MULTICAST,
        PACKET_OTHERHOST as PACKET_OTHERHOST,
        PACKET_OUTGOING as PACKET_OUTGOING,
        PF_CAN as PF_CAN,
        PF_PACKET as PF_PACKET,
        PF_RDS as PF_RDS,
        RDS_CANCEL_SENT_TO as RDS_CANCEL_SENT_TO,
        RDS_CMSG_RDMA_ARGS as RDS_CMSG_RDMA_ARGS,
        RDS_CMSG_RDMA_DEST as RDS_CMSG_RDMA_DEST,
        RDS_CMSG_RDMA_MAP as RDS_CMSG_RDMA_MAP,
        RDS_CMSG_RDMA_STATUS as RDS_CMSG_RDMA_STATUS,
        RDS_CMSG_RDMA_UPDATE as RDS_CMSG_RDMA_UPDATE,
        RDS_CONG_MONITOR as RDS_CONG_MONITOR,
        RDS_FREE_MR as RDS_FREE_MR,
        RDS_GET_MR as RDS_GET_MR,
        RDS_GET_MR_FOR_DEST as RDS_GET_MR_FOR_DEST,
        RDS_RDMA_DONTWAIT as RDS_RDMA_DONTWAIT,
        RDS_RDMA_FENCE as RDS_RDMA_FENCE,
        RDS_RDMA_INVALIDATE as RDS_RDMA_INVALIDATE,
        RDS_RDMA_NOTIFY_ME as RDS_RDMA_NOTIFY_ME,
        RDS_RDMA_READWRITE as RDS_RDMA_READWRITE,
        RDS_RDMA_SILENT as RDS_RDMA_SILENT,
        RDS_RDMA_USE_ONCE as RDS_RDMA_USE_ONCE,
        RDS_RECVERR as RDS_RECVERR,
        SOL_ALG as SOL_ALG,
        SOL_CAN_BASE as SOL_CAN_BASE,
        SOL_CAN_RAW as SOL_CAN_RAW,
        SOL_RDS as SOL_RDS,
        SOL_TIPC as SOL_TIPC,
        TIPC_ADDR_ID as TIPC_ADDR_ID,
        TIPC_ADDR_NAME as TIPC_ADDR_NAME,
        TIPC_ADDR_NAMESEQ as TIPC_ADDR_NAMESEQ,
        TIPC_CFG_SRV as TIPC_CFG_SRV,
        TIPC_CLUSTER_SCOPE as TIPC_CLUSTER_SCOPE,
        TIPC_CONN_TIMEOUT as TIPC_CONN_TIMEOUT,
        TIPC_CRITICAL_IMPORTANCE as TIPC_CRITICAL_IMPORTANCE,
        TIPC_DEST_DROPPABLE as TIPC_DEST_DROPPABLE,
        TIPC_HIGH_IMPORTANCE as TIPC_HIGH_IMPORTANCE,
        TIPC_IMPORTANCE as TIPC_IMPORTANCE,
        TIPC_LOW_IMPORTANCE as TIPC_LOW_IMPORTANCE,
        TIPC_MEDIUM_IMPORTANCE as TIPC_MEDIUM_IMPORTANCE,
        TIPC_NODE_SCOPE as TIPC_NODE_SCOPE,
        TIPC_PUBLISHED as TIPC_PUBLISHED,
        TIPC_SRC_DROPPABLE as TIPC_SRC_DROPPABLE,
        TIPC_SUB_CANCEL as TIPC_SUB_CANCEL,
        TIPC_SUB_PORTS as TIPC_SUB_PORTS,
        TIPC_SUB_SERVICE as TIPC_SUB_SERVICE,
        TIPC_SUBSCR_TIMEOUT as TIPC_SUBSCR_TIMEOUT,
        TIPC_TOP_SRV as TIPC_TOP_SRV,
        TIPC_WAIT_FOREVER as TIPC_WAIT_FOREVER,
        TIPC_WITHDRAWN as TIPC_WITHDRAWN,
        TIPC_ZONE_SCOPE as TIPC_ZONE_SCOPE,
    )
if sys.platform == "linux" and sys.version_info >= (3, 7):
    from _socket import (
        CAN_ISOTP as CAN_ISOTP,
        IOCTL_VM_SOCKETS_GET_LOCAL_CID as IOCTL_VM_SOCKETS_GET_LOCAL_CID,
        SO_VM_SOCKETS_BUFFER_MAX_SIZE as SO_VM_SOCKETS_BUFFER_MAX_SIZE,
        SO_VM_SOCKETS_BUFFER_MIN_SIZE as SO_VM_SOCKETS_BUFFER_MIN_SIZE,
        SO_VM_SOCKETS_BUFFER_SIZE as SO_VM_SOCKETS_BUFFER_SIZE,
        TCP_NOTSENT_LOWAT as TCP_NOTSENT_LOWAT,
        VM_SOCKETS_INVALID_VERSION as VM_SOCKETS_INVALID_VERSION,
        VMADDR_CID_ANY as VMADDR_CID_ANY,
        VMADDR_CID_HOST as VMADDR_CID_HOST,
        VMADDR_PORT_ANY as VMADDR_PORT_ANY,
    )
if sys.platform == "linux" and sys.version_info >= (3, 8):
    from _socket import (
        CAN_BCM_CAN_FD_FRAME as CAN_BCM_CAN_FD_FRAME,
        CAN_BCM_RX_ANNOUNCE_RESUME as CAN_BCM_RX_ANNOUNCE_RESUME,
        CAN_BCM_RX_CHECK_DLC as CAN_BCM_RX_CHECK_DLC,
        CAN_BCM_RX_FILTER_ID as CAN_BCM_RX_FILTER_ID,
        CAN_BCM_RX_NO_AUTOTIMER as CAN_BCM_RX_NO_AUTOTIMER,
        CAN_BCM_RX_RTR_FRAME as CAN_BCM_RX_RTR_FRAME,
        CAN_BCM_SETTIMER as CAN_BCM_SETTIMER,
        CAN_BCM_STARTTIMER as CAN_BCM_STARTTIMER,
        CAN_BCM_TX_ANNOUNCE as CAN_BCM_TX_ANNOUNCE,
        CAN_BCM_TX_COUNTEVT as CAN_BCM_TX_COUNTEVT,
        CAN_BCM_TX_CP_CAN_ID as CAN_BCM_TX_CP_CAN_ID,
        CAN_BCM_TX_RESET_MULTI_IDX as CAN_BCM_TX_RESET_MULTI_IDX,
    )
if sys.platform == "linux" and sys.version_info >= (3, 9):
    from _socket import (
        CAN_J1939 as CAN_J1939,
        CAN_RAW_JOIN_FILTERS as CAN_RAW_JOIN_FILTERS,
        J1939_EE_INFO_NONE as J1939_EE_INFO_NONE,
        J1939_EE_INFO_TX_ABORT as J1939_EE_INFO_TX_ABORT,
        J1939_FILTER_MAX as J1939_FILTER_MAX,
        J1939_IDLE_ADDR as J1939_IDLE_ADDR,
        J1939_MAX_UNICAST_ADDR as J1939_MAX_UNICAST_ADDR,
        J1939_NLA_BYTES_ACKED as J1939_NLA_BYTES_ACKED,
        J1939_NLA_PAD as J1939_NLA_PAD,
        J1939_NO_ADDR as J1939_NO_ADDR,
        J1939_NO_NAME as J1939_NO_NAME,
        J1939_NO_PGN as J1939_NO_PGN,
        J1939_PGN_ADDRESS_CLAIMED as J1939_PGN_ADDRESS_CLAIMED,
        J1939_PGN_ADDRESS_COMMANDED as J1939_PGN_ADDRESS_COMMANDED,
        J1939_PGN_MAX as J1939_PGN_MAX,
        J1939_PGN_PDU1_MAX as J1939_PGN_PDU1_MAX,
        J1939_PGN_REQUEST as J1939_PGN_REQUEST,
        SCM_J1939_DEST_ADDR as SCM_J1939_DEST_ADDR,
        SCM_J1939_DEST_NAME as SCM_J1939_DEST_NAME,
        SCM_J1939_ERRQUEUE as SCM_J1939_ERRQUEUE,
        SCM_J1939_PRIO as SCM_J1939_PRIO,
        SO_J1939_ERRQUEUE as SO_J1939_ERRQUEUE,
        SO_J1939_FILTER as SO_J1939_FILTER,
        SO_J1939_PROMISC as SO_J1939_PROMISC,
        SO_J1939_SEND_PRIO as SO_J1939_SEND_PRIO,
    )
if sys.platform == "win32":
    from _socket import (
        RCVALL_IPLEVEL as RCVALL_IPLEVEL,
        RCVALL_MAX as RCVALL_MAX,
        RCVALL_OFF as RCVALL_OFF,
        RCVALL_ON as RCVALL_ON,
        RCVALL_SOCKETLEVELONLY as RCVALL_SOCKETLEVELONLY,
        SIO_KEEPALIVE_VALS as SIO_KEEPALIVE_VALS,
        SIO_LOOPBACK_FAST_PATH as SIO_LOOPBACK_FAST_PATH,
        SIO_RCVALL as SIO_RCVALL,
    )

_T = TypeVar("_T")

# Re-exported from errno
EBADF: int
EAGAIN: int
EWOULDBLOCK: int

class AddressFamily(IntEnum):
    AF_UNIX: int
    AF_INET: int
    AF_INET6: int
    AF_AAL5: int
    AF_ALG: int
    AF_APPLETALK: int
    AF_ASH: int
    AF_ATMPVC: int
    AF_ATMSVC: int
    AF_AX25: int
    AF_BLUETOOTH: int
    AF_BRIDGE: int
    AF_CAN: int
    AF_DECnet: int
    AF_ECONET: int
    AF_IPX: int
    AF_IRDA: int
    AF_KEY: int
    AF_LINK: int
    AF_LLC: int
    AF_NETBEUI: int
    AF_NETLINK: int
    AF_NETROM: int
    AF_PACKET: int
    AF_PPPOX: int
    AF_QIPCRTR: int
    AF_RDS: int
    AF_ROSE: int
    AF_ROUTE: int
    AF_SECURITY: int
    AF_SNA: int
    AF_SYSTEM: int
    AF_TIPC: int
    AF_UNSPEC: int
    AF_VSOCK: int
    AF_WANPIPE: int
    AF_X25: int

AF_UNIX: AddressFamily
AF_INET: AddressFamily
AF_INET6: AddressFamily
AF_AAL5: AddressFamily
AF_APPLETALK: AddressFamily
AF_ASH: AddressFamily
AF_ATMPVC: AddressFamily
AF_ATMSVC: AddressFamily
AF_AX25: AddressFamily
AF_BRIDGE: AddressFamily
AF_DECnet: AddressFamily
AF_ECONET: AddressFamily
AF_IPX: AddressFamily
AF_IRDA: AddressFamily
AF_KEY: AddressFamily
AF_LLC: AddressFamily
AF_NETBEUI: AddressFamily
AF_NETROM: AddressFamily
AF_PPPOX: AddressFamily
AF_ROSE: AddressFamily
AF_ROUTE: AddressFamily
AF_SECURITY: AddressFamily
AF_SNA: AddressFamily
AF_SYSTEM: AddressFamily
AF_UNSPEC: AddressFamily
AF_WANPIPE: AddressFamily
AF_X25: AddressFamily
if sys.platform == "linux":
    AF_CAN: AddressFamily
    AF_PACKET: AddressFamily
    AF_RDS: AddressFamily
    AF_TIPC: AddressFamily
    AF_ALG: AddressFamily
    AF_NETLINK: AddressFamily
    if sys.version_info >= (3, 7):
        AF_VSOCK: AddressFamily
    if sys.version_info >= (3, 8):
        AF_QIPCRTR: AddressFamily
AF_LINK: AddressFamily  # availability: BSD, macOS
if sys.platform != "win32" and sys.platform != "darwin":
    AF_BLUETOOTH: AddressFamily

class SocketKind(IntEnum):
    SOCK_STREAM: int
    SOCK_DGRAM: int
    SOCK_RAW: int
    SOCK_RDM: int
    SOCK_SEQPACKET: int
    SOCK_CLOEXEC: int
    SOCK_NONBLOCK: int

SOCK_STREAM: SocketKind
SOCK_DGRAM: SocketKind
SOCK_RAW: SocketKind
SOCK_RDM: SocketKind
SOCK_SEQPACKET: SocketKind
if sys.platform == "linux":
    SOCK_CLOEXEC: SocketKind
    SOCK_NONBLOCK: SocketKind

class MsgFlag(IntFlag):
    MSG_CTRUNC: int
    MSG_DONTROUTE: int
    MSG_DONTWAIT: int
    MSG_EOR: int
    MSG_OOB: int
    MSG_PEEK: int
    MSG_TRUNC: int
    MSG_WAITALL: int

MSG_BCAST: MsgFlag
MSG_BTAG: MsgFlag
MSG_CMSG_CLOEXEC: MsgFlag
MSG_CONFIRM: MsgFlag
MSG_CTRUNC: MsgFlag
MSG_DONTROUTE: MsgFlag
MSG_DONTWAIT: MsgFlag
MSG_EOF: MsgFlag
MSG_EOR: MsgFlag
MSG_ERRQUEUE: MsgFlag
MSG_ETAG: MsgFlag
MSG_FASTOPEN: MsgFlag
MSG_MCAST: MsgFlag
MSG_MORE: MsgFlag
MSG_NOSIGNAL: MsgFlag
MSG_NOTIFICATION: MsgFlag
MSG_OOB: MsgFlag
MSG_PEEK: MsgFlag
MSG_TRUNC: MsgFlag
MSG_WAITALL: MsgFlag

class AddressInfo(IntFlag):
    AI_ADDRCONFIG: int
    AI_ALL: int
    AI_CANONNAME: int
    AI_NUMERICHOST: int
    AI_NUMERICSERV: int
    AI_PASSIVE: int
    AI_V4MAPPED: int

AI_ADDRCONFIG: AddressInfo
AI_ALL: AddressInfo
AI_CANONNAME: AddressInfo
AI_DEFAULT: AddressInfo
AI_MASK: AddressInfo
AI_NUMERICHOST: AddressInfo
AI_NUMERICSERV: AddressInfo
AI_PASSIVE: AddressInfo
AI_V4MAPPED: AddressInfo
AI_V4MAPPED_CFG: AddressInfo

if sys.platform == "win32":
    errorTab: dict[int, str]  # undocumented

class socket(_socket.socket):
    def __init__(
        self, family: AddressFamily | int = ..., type: SocketKind | int = ..., proto: int = ..., fileno: int | None = ...
    ) -> None: ...
    def __enter__(self: Self) -> Self: ...
    def __exit__(self, *args: object) -> None: ...
    def dup(self: _T) -> _T: ...  # noqa: F811
    def accept(self) -> tuple[socket, _RetAddress]: ...
    # Note that the makefile's documented windows-specific behavior is not represented
    # mode strings with duplicates are intentionally excluded
    @overload
    def makefile(
        self,
        mode: Literal["r", "w", "rw", "wr", ""] = ...,
        buffering: int | None = ...,
        *,
        encoding: str | None = ...,
        errors: str | None = ...,
        newline: str | None = ...,
    ) -> TextIO: ...
    @overload
    def makefile(
        self,
        mode: Literal["b", "rb", "br", "wb", "bw", "rwb", "rbw", "wrb", "wbr", "brw", "bwr"],
        buffering: int | None = ...,
        *,
        encoding: str | None = ...,
        errors: str | None = ...,
        newline: str | None = ...,
    ) -> BinaryIO: ...
    def sendfile(self, file: BinaryIO, offset: int = ..., count: int | None = ...) -> int: ...
    @property
    def family(self) -> AddressFamily: ...  # type: ignore
    @property
    def type(self) -> SocketKind: ...  # type: ignore
    def get_inheritable(self) -> bool: ...
    def set_inheritable(self, inheritable: bool) -> None: ...

def fromfd(fd: _FD, family: AddressFamily | int, type: SocketKind | int, proto: int = ...) -> socket: ...

if sys.platform != "win32":
    if sys.version_info >= (3, 9):
        # flags and address appear to be unused in send_fds and recv_fds
        def send_fds(
            sock: socket, buffers: Iterable[bytes], fds: bytes | Iterable[int], flags: int = ..., address: None = ...
        ) -> int: ...
        def recv_fds(sock: socket, bufsize: int, maxfds: int, flags: int = ...) -> tuple[bytes, list[int], int, Any]: ...

if sys.platform == "win32":
    def fromshare(info: bytes) -> socket: ...

if sys.platform == "win32":
    def socketpair(family: int = ..., type: int = ..., proto: int = ...) -> tuple[socket, socket]: ...

else:
    def socketpair(  # type: ignore
        family: int | AddressFamily | None = ..., type: SocketType | int = ..., proto: int = ...
    ) -> tuple[socket, socket]: ...

class SocketIO(RawIOBase):
    def __init__(self, sock: socket, mode: Literal["r", "w", "rw", "rb", "wb", "rwb"]) -> None: ...
    def readinto(self, b: WriteableBuffer) -> int | None: ...
    def write(self, b: ReadableBuffer) -> int | None: ...
    @property
    def name(self) -> int: ...  # return value is really "int"
    @property
    def mode(self) -> Literal["rb", "wb", "rwb"]: ...

def getfqdn(name: str = ...) -> str: ...
def create_connection(
    address: tuple[str | None, int],
    timeout: float | None = ...,  # noqa: F811
    source_address: tuple[bytearray | bytes | str, int] | None = ...,
) -> socket: ...

if sys.version_info >= (3, 8):
    def has_dualstack_ipv6() -> bool: ...
    def create_server(
        address: _Address, *, family: int = ..., backlog: int | None = ..., reuse_port: bool = ..., dualstack_ipv6: bool = ...
    ) -> socket: ...

# the 5th tuple item is an address
def getaddrinfo(
    host: bytearray | bytes | str | None,
    port: str | int | None,
    family: int = ...,
    type: int = ...,
    proto: int = ...,
    flags: int = ...,
) -> list[tuple[AddressFamily, SocketKind, int, str, tuple[str, int] | tuple[str, int, int, int]]]: ...
